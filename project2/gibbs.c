#include <math.h>
#include <time.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>


void Gibbs(const double **filters, const int *h, const int *w, const int H, const int W, const int num_filters, const int num_bins, const double **hists_gt, const double *bounds, const int sweep, const int is_frame, const double **multipliers, int *img_syn, double **responses_syn, double **hists_syn)
{
	time_t start;
	time_t end;

	int s, i, j, k, l, b;
	int h_filter, w_filter, x, y;
	int gray_ori, gray_new;
	int modified_response_x, modified_response_y;
	int minimum;

	double T;
	double lower_bound, upper_bound, bin_width;
	
	double response_z1, response_z2;
	int response_bin1, response_bin2;

	double *probs = (double *)malloc(8 * sizeof(double));
	double ***hists_tmp = (double ***)malloc(8 * sizeof(double **));
	for (k = 0; k < 8; k++)
	{
		hists_tmp[k] = (double **)malloc(num_filters * sizeof(double *));
		for (l = 0; l < num_filters; l++)
			hists_tmp[k][l] = (double *)malloc(num_bins * sizeof(double));
	}
	
	double *band_weights = (double *)malloc(num_bins * sizeof(double));
	for (b = 0; b < num_bins; b++)
		band_weights[b] = (double)(1 + fabs((num_bins-1)/2 - b));

	double sum_probs;
	double random_num;
	double sum_error, mean_error, sum_hist_syn, sum_hist_gt;

	mean_error = 0;
	for (l = 0; l < num_filters; l++)
	{
		sum_error = 0;
		sum_hist_syn = 0;
		sum_hist_gt = 0;
		printf("hist_syn: ");
		for (b = 0; b < num_bins; b++)
		{
			printf("%f ", hists_syn[l][b]);
			sum_hist_syn += hists_syn[l][b];
			sum_error += fabs(hists_syn[l][b]-hists_gt[l][b]);
		}
		printf("= %f\nhist_gt : ", sum_hist_syn);
		for (b = 0; b < num_bins; b++)
		{
			sum_hist_gt += hists_gt[l][b];
			printf("%f ", hists_gt[l][b]);
		}
		printf("= %f\n\n", sum_hist_gt);
		mean_error += sum_error;
	}
	mean_error /= num_filters;
	printf("Gibbs iteration 0, error = %.6lf\n\n\n", mean_error);

	T = 0.1;
	s = 0;
    
    while (mean_error > 0.01 && s < 50)
    {
        time(&start);
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++)
			{
				gray_ori = img_syn[W*i+j];

				for (k = 0; k < 8; k++)
					probs[k] = 0;				
				// compute conditional probabilities
				for (k = 0; k < 8; k++)
					for (l = 0; l < num_filters; l++)
					{
						for (b = 0; b < num_bins; b++)
							hists_tmp[k][l][b] = hists_syn[l][b] * H * W;

						lower_bound = bounds[2*l];
						upper_bound = bounds[2*l+1];
						bin_width = (upper_bound - lower_bound) / num_bins;
						h_filter = h[l];
						w_filter = w[l];
						for (y = 0; y < h_filter; y++)
							for (x = 0; x < w_filter; x++)
							{
								// `current_pos` multiplied by pixel (y, x) on filter
								modified_response_x = j + (w_filter-1)/2 - x;
								modified_response_y = i + (h_filter-1)/2 - y;
								
								// in case of boundary
								if (modified_response_x < 0 || modified_response_x > W-1 || modified_response_y < 0 || modified_response_y > H-1)
									continue;
								
								response_z1 = responses_syn[l][W * modified_response_y + modified_response_x];
								response_bin1 = (int)(floor((response_z1-lower_bound) / bin_width));
								if (response_bin1 < 0)
									response_bin1 = 0;
								else if (response_bin1 > num_bins - 1)
									response_bin1 = num_bins - 1;
								hists_tmp[k][l][response_bin1] -= 1.0;

								response_z2 = response_z1 + (double)( (k-gray_ori) * filters[l][w_filter*y+x] );
								response_bin2 = (int)(floor((response_z2-lower_bound) / bin_width));
								if (response_bin2 < 0)
									response_bin2 = 0;
								else if (response_bin2 > num_bins - 1)
									response_bin2 = num_bins - 1;
								hists_tmp[k][l][response_bin2] += 1.0;
								// printf("lower_bound=%f: (%f, %d) + (%d - %d) * %f (bin_width=%f) -> (%f, %d)\n", lower_bound, response_z1, response_bin1, k, gray_ori, filters[l][w_filter*y+x], bin_width, response_z2, response_bin2);
							}

						for (b = 0; b < num_bins; b++)
						{
							if (is_frame == 0)
								// Julesz ensemble
								probs[k] += ( fabs( hists_tmp[k][l][b] - hists_gt[l][b]*H*W ) * band_weights[b] );
							else
								// FRAME
								probs[k] += ( hists_tmp[k][l][b] * multipliers[l][b] );
							
							hists_tmp[k][l][b] /= (double)(H * W);
						}
					}

				minimum = INT_MAX;
				for (k = 0; k < 8; k++)
					if (probs[k] < minimum)
						minimum = probs[k];

				sum_probs = 0;
				for (k = 0; k < 8; k++)
				{
					probs[k] = probs[k] - minimum;
					sum_probs += probs[k];
				}

				for (k = 0; k < 8; k++)
				{
					probs[k] = probs[k] / (sum_probs+1e-6);
					probs[k] = exp(-probs[k] / T);
				}

				sum_probs = 0;
				for (k = 0; k < 8; k++)
					sum_probs += probs[k];

				for (k = 0; k < 8; k++)
					probs[k] /= sum_probs;
				
				// sample
				random_num = (double)rand() / RAND_MAX;
				for (k = 0; k < 8; k++)
				{
					if (random_num < probs[k])
					{
						gray_new = k;
						break;
					}
					else
						random_num = random_num - probs[k];
				}

				img_syn[W*i+j] = gray_new;

				// update
				for (l = 0; l < num_filters; l++)
				{
					lower_bound = bounds[2*l];
					upper_bound = bounds[2*l+1];
					bin_width = (upper_bound - lower_bound) / num_bins;
					h_filter = h[l];
					w_filter = w[l];
					
					// update responses
					for (y = 0; y < h_filter; y++)
						for (x = 0; x < w_filter; x++)
						{
							// `current_pos` multiplied by pixel (y, x) on filter
							modified_response_x = j + (int)((w_filter-1)/2) - x;
							modified_response_y = i + (int)((h_filter-1)/2) - y;
							
							// in case of boundary
							if (modified_response_x < 0 || modified_response_x > W-1 || modified_response_y < 0 || modified_response_y > H-1)
								continue;

							responses_syn[l][W*modified_response_y+modified_response_x] += (double)( (gray_new-gray_ori) * filters[l][w_filter*y+x] );
						}
					
					// update histograms
					for (b = 0; b < num_bins; b++)
						hists_syn[l][b] = hists_tmp[gray_new][l][b];
				}
			}
        
        T *= 0.96;
        s += 1;

        mean_error = 0;
		for (l = 0; l < num_filters; l++)
		{
			sum_error = 0;
			for (b = 0; b < num_bins; b++)
				sum_error += fabs(hists_syn[l][b]-hists_gt[l][b]);
			mean_error += sum_error;
		}
		mean_error /= num_filters;

        time(&end);
        if (s % 10 == 0)
		{
            printf("Gibbs iteration %d took %.2lf s, error = %.6lf\n", s, difftime(end, start), mean_error);
			for (l = 0; l < num_filters; l++)
			{
				sum_hist_syn = 0;
				sum_hist_gt = 0;
				printf("hist_syn: ");
				for (b = 0; b < num_bins; b++)
				{
					printf("%f ", hists_syn[l][b]);
					sum_hist_syn += hists_syn[l][b];
				}
				printf("= %f\nhist_gt : ", sum_hist_syn);
				for (b = 0; b < num_bins; b++)
				{
					sum_hist_gt += hists_gt[l][b];
					printf("%f ", hists_gt[l][b]);
				}
				printf("= %f\n\n\n", sum_hist_gt);
			}
		}
    }

	// release memory

	free(probs);
	free(band_weights);

	for (k = 0; k < 8; k++)
	{
		for (l = 0; l < num_filters; l++)
			free(hists_tmp[k][l]);
		free(hists_tmp[k]);
	}
	free(hists_tmp);
}


void test_int_1D(int *input, int len)
{
	printf("1D int in C: ");
	for (int i = 0; i < len; i++)
		printf("%d ", input[i]);
	printf("\n");
}


void test_double_1D(double *input, int len)
{
	printf("1D double in C: ");
	for (int i = 0; i < len; i++)
		printf("%f ", input[i]);
	printf("\n");
}


void test_double_2D(double **input, int len1, int len2)
{
	printf("2D double in C: ");
	for (int i = 0; i < len1; i++)
	{
		for (int j = 0; j < len2; j++)
			printf("%f ", input[i][j]);
		printf("\n");
	}
}

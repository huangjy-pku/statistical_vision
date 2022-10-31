all: lib_gibbs.so

lib_gibbs.so: gibbs.c
	$(CC) -Wall -g -fPIC -shared -o $@ $? -lc

clean:
	rm -f test lib_gibbs.o lib_gibbs.so *.pyc
	rm -rf build
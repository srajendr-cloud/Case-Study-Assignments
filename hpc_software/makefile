CC=mpicc
CFLAGS=-O2

all: q1_sendrecv q1_collective q2_test q3_norm

q1_sendrecv:
	$(CC) $(CFLAGS) q1_sendrecv.c -o q1_sendrecv

q1_collective:
	$(CC) $(CFLAGS) q1_collective.c -o q1_collective

q2_test:
	$(CC) $(CFLAGS) q2_test.c decomp1d.c -o q2_test

q3_norm:
	$(CC) $(CFLAGS) q3_norm.c decomp1d.c -o q3_norm

clean:
	rm -f q1_sendrecv q1_collective q2_test q3_norm


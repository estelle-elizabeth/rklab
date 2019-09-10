#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/***********************************************************
 Implementation of bloom filter goes here 
 **********************************************************/

#include "bloom.h"

/* Constants for bloom filter implementation */
const int H1PRIME = 4189793;
const int H2PRIME = 3296731;
const int BLOOM_HASH_NUM = 10;

long long BIG_PRIME = 5003943032159437;
long long int asc = 256;

const int PRINT_RK_HASH = 5;
const int PRINT_BLOOM_BITS = 160;

long long
madd(long long a, long long b)
{
	return ((a+b)>BIG_PRIME?(a+b-BIG_PRIME):(a+b));
}

/* modulo substraction */
long long
mdel(long long a, long long b)
{
	return ((a>b)?(a-b):(a+BIG_PRIME-b));
}

/* modulo multiplication*/
long long
mmul(long long a, long long b)
{
	return ((a*b) % BIG_PRIME);
}

/* The hash function used by the bloom filter */
int
hash_i(int i, /* which of the BLOOM_HASH_NUM hashes to use */ 
       long long x /* a long long value to be hashed */)
{
	return ((x % H1PRIME) + i*(x % H2PRIME) + 1 + i*i);
}

/* Initialize a bloom filter by allocating a character array that can pack bsz bits.
   (each char represents 8 bits)
   Furthermore, clear all bits for the allocated character array. 
   Hint:  use the malloc and bzero library function 
	 Return value is the newly initialized bloom_filter struct.*/
bloom_filter 
bloom_init(int bsz /* size of bitmap to allocate in bits*/ )
{
	bloom_filter f;
	f.bsz = bsz;

	/* your code here*/
	int size = bsz/8;
	if ((bsz % 8) > 0) size++;

	f.buf = (char *) malloc(size);
	bzero(f.buf, size);

	// char bit = 1;
	// for (int i = 0; i < size*8; i++){
	// 	int div = i/8;
	// 	int rem = 7-i%8;
	// 	// char *byte = &f.buf[div];
	// 	bit = bit << rem;
	// 	printf("%d\n", i);
	// 	printf("%d\n", bit & f.buf[div]);
	// }
	// memset(f.buf, 0, size);
	return f;
}

/* Add elm into the given bloom filter*/
void
bloom_add(bloom_filter f,
          long long elm /* the element to be added (a RK hash value) */)
{
	int hashed = 0, div = 0, rem = 0;
	char bit;
	int size = f.bsz;
	for (int i = 0; i < BLOOM_HASH_NUM; i++){
		hashed = hash_i(i, elm) % size;
		div = hashed/8;
		rem = 7-(hashed%8); //for big_endian
		// printf("hash is %d\n", div);
		// char *byte = &f.buf[div];
		bit = 1 << rem;

		// printf("%s\n", "starting");
		// char a = 1;
		// for (int i = 0; i < sizeof(int)*8; i++){
		// 	// int c = i/8;
		// 	// int d = 7-(i%8);
		char *byte = &f.buf[div];
		// 	a = a<< rem;
		// 	// printf("%d\n", i);
		// 	// printf("%d\n", a & bit);
		// }

		*byte = *byte | bit;
		// printf("%s\n", &f.buf[div]);
		// f.buf[div] = f.buf[div] | bit;

	}
	// printf("%s\n", "starting");
	// char a = 1;
	// for (int i = 0; i < size; i++){
	// 	int c = i/8;
	// 	int d = 7-(i%8);
	// 	// char *byte = &f.buf[div];
	// 	a = a<< d;
	// 	// printf("%d\n", i);
	// 	printf("%d\n", a & f.buf[c]);
	// }

	return;
}

/* Query if elm is probably in the given bloom filter */ 
int
bloom_query(bloom_filter f,
            long long elm /* the query element */ )
{	
	int hashed = 0, div = 0, rem = 0;
	char bit;
	int count = 0, size = f.bsz;
	for (int i = 0; i < BLOOM_HASH_NUM; i++){
		hashed = hash_i(i, elm)%size;
		div = hashed/8;
		rem = 7-(hashed%8);
		// char *byte = &f.buf[div];
		bit = 1 << rem;
		// if (bit & *byte) count++;
		if (bit & f.buf[div]) count++;
		else{
			break;
		}
	}

	if (count == BLOOM_HASH_NUM) return 1;

	return 0;
}


void 
bloom_free(bloom_filter *f)
{
	free(f->buf);
	f->buf = f->bsz = 0;
}

/* print out the first count bits in the bloom filter */
void
bloom_print(bloom_filter f,
            int count     /* number of bits to display*/ )
{
	int i;

	assert(count % 8 == 0);

	for(i=0; i< (f.bsz>>3) && i < (count>>3); i++) {
		printf("%02x ", (unsigned char)(f.buf[i]));
	}
	printf("\n");
	return;
}


int
rabin_karp_batchmatch(int bsz,        /* size of bitmap (in bits) to be used */
                      int k,          /* chunk length to be matched */
                      const char *qs, /* query docoument (X)*/
                      int m,          /* query document length */ 
                      const char *ts, /* to-be-matched document (Y) */
                      int n           /* to-be-matched document length*/)
{
	bloom_filter char_array = bloom_init(bsz);


	long long int hashqs = 0, hashts = 0;
	long long int pow, largest = 0;
	int chunk = m/k, count = 0;
	for (int i = 0; i < m; i+=k){
		
		
		for (int a = 0; a < k; a++){
			pow = 1;
			for (int j = 1; j < k-a; j++){
				pow = mmul(pow, asc);
			}

			if (a==0) {
				largest = pow;
			}
			hashqs = madd(hashqs, mmul(pow, qs[i]));
			if (i==0){
				hashts = madd(hashts, mmul(pow, ts[i]));
			}
			
		}
		bloom_add(char_array, hashqs);
	}

	bloom_print(char_array, PRINT_BLOOM_BITS);
	// printf("\n");

	// printf("count is%d\n", count);
	for (int i = 0; (i+k) <= n; i++){
			if (i > 0){
				hashts = mdel(hashts, mmul(largest, ts[i-1]));
				hashts = mmul(hashts, asc);
				hashts = madd(hashts, ts[i+k-1]);
			}

			// printf("bloom query is %d\n", bloom_query(char_array, hashts));
			// printf("i s%d\n", i+k);
			if (bloom_query(char_array, hashts)){
				// for(int index=0; index<m; index+=k){
				// 	if(strncmp((qs+index), ts+i, k)==0){
				// 		printf("%s\n", "count is added");
				// 	  	count +=1;
				// 	}
				// }
				// // count++;
				int match;
				for(int j=0; j<m; j+=k){
					match = 0;
					for (int a = 0; a<k; a++){
						if (qs[a] == ts[a+i]){
							match++;
						}

						else{
							break;
						}

					}

					if (match == k) count++;
					
				}
			}

			
	
	
}
	// printf("count is%d\n", count);
	bloom_free(&char_array);
	return count;
}


int main(){
	int bsz = 200;
	char qs[] = "abdc";
	char ts[] = "dabdc";
	int count = (rabin_karp_batchmatch(bsz, 2, qs, strlen(qs), ts, strlen(ts)));
	// printf("count is %d\n", count);

	int size = 10;
	char *buf = (char *) malloc(size);
	bzero(buf, size);

	int hashed = 0, div = 0, rem = 0;
	char bit;
	bit = 1 << 3;
	char *byte = &buf[4];
	buf[4] = buf[4] | bit;
	printf("read %d\n", (1<<3) & buf[4]);
	free(buf);
	// for (int i = 0; i < BLOOM_HASH_NUM; i++){
	// 	hashed = hash_i(i, elm) % size;
	// 	div = hashed/8;
	// 	rem = 7-(hashed%8); //for big_endian
	// 	// printf("hash is %d\n", div);
	// 	// char *byte = &f.buf[div];
	// 	bit = 1 << rem;

	// 	// printf("%s\n", "starting");
	// 	// char a = 1;
	// 	// for (int i = 0; i < sizeof(int)*8; i++){
	// 	// 	// int c = i/8;
	// 	// 	// int d = 7-(i%8);
	// 	char *byte = &f.buf[div];
	// 	// 	a = a<< rem;
	// 	// 	// printf("%d\n", i);
	// 	// 	// printf("%d\n", a & bit);
	// 	// }

	// 	*byte = *byte | bit;
	// 	// printf("%s\n", &f.buf[div]);
	// 	// f.buf[div] = f.buf[div] | bit;

	// }

	// for (int i = 0; (i+2) < 4; i++){
	// 	printf("%s\n", "here");
	// }

	// int bit = 1;
	// for (int i = 0; i < BLOOM_HASH_NUM; i++){
	// 	bit = 15;

		// char a = 1;
		// for (int i = 0; i < sizeof(int)*8; i++){
		// 	// int c = i/8;
		// 	// int d = 7-(i%8);
		// 	// char *byte = &f.buf[div];
		// 	a = a<< 3;
		// 	// printf("%d\n", i);
		// 	printf("%d\n", a & bit);
		// }

		// if (15 & (1<<0)){
		// 	printf("%s\n", "yes");
		// }
		// else printf("%s\n", "no");

		// int num = 7;
		// num = num | (1<<2);
		// printf("%d\n", num);
	// }

	return 0;
}




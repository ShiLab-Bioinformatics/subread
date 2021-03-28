/*--------------------------------------------------------------------------*\
 *				   -----===== HashTable =====-----
 *
 * Author: Keith Pomakis (pomakis@pobox.com)
 * Date:   August, 1998
 * Released to the public domain.
 *
 *--------------------------------------------------------------------------
 * $Id: hashtable.c,v 9999.42 2019/12/05 03:05:37 cvs Exp $
\*--------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include "hashtable.h"
#include "core.h"
#include "input-files.h"

static int pointercmp(const void *pointer1, const void *pointer2);
static srUInt_64 pointerHashFunction(const void *pointer);
static int isProbablePrime(srInt_64 number);
static srInt_64 calculateIdealNumOfBuckets(HashTable *hashTable);

srInt_64 long_random_val(){
	srInt_64 ret = 0;
	if(RAND_MAX<255){
		SUBREADprintf("Is this a embedded computer????\n");
		return -1;
	}
	int i;
	for(i=0;i<8;i++){
		if(i>0)ret = (ret << 8) ^ (myrand_rand() & 0xff);
		else ret = (ret << 8) ^ (myrand_rand() & 0x7f);
	}
	return ret;
}

HashTable * ArrayListToLookupTable_Int(ArrayList * arr){
	srInt_64 x1;
	HashTable * ret = HashTableCreate( max(1,arr -> numOfElements / 6) );
	for(x1 = 0; x1 < arr->numOfElements; x1++) HashTablePut(ret, ArrayListGet(arr, x1)+1, NULL+x1+1);
	return ret;
}

int ArrayListContainsPtr(ArrayList * list, void * who){
	srInt_64 x1;
	for(x1 = 0; x1 < list->numOfElements; x1++) if(list->elementList [ x1 ]==who) return 1;
	return 0;
}

void * ArrayListShift(ArrayList * list){
	if(list->numOfElements<1) return NULL;
	void *ret = list->elementList [0];
	srInt_64 xx;
	list->numOfElements -- ;
	for(xx=0; xx<list->numOfElements; xx++) list->elementList [ xx ] = list->elementList [ xx+1 ];
	return ret;
}

void * ArrayListPop(ArrayList * list){
	if(list->numOfElements<1) return NULL;
	return list->elementList [ -- list->numOfElements];
}


void * ArrayListRandom(ArrayList * list){
	srInt_64 ii = long_random_val() % list -> numOfElements;
	return list -> elementList[ii];
}


ArrayList * ArrayListCreate(int init_capacity){
	ArrayList * ret = malloc(sizeof(ArrayList));
	memset(ret,0,sizeof(ArrayList));
	ret -> capacityOfElements = init_capacity;
	ret -> elementList = malloc(sizeof(void *)*init_capacity);
	return ret;
}

ArrayList * ArrayListDuplicate(ArrayList * ori){
	ArrayList * ret = malloc(sizeof(ArrayList));
	memset(ret,0,sizeof(ArrayList));
	ret -> capacityOfElements = ori -> numOfElements;
	ret -> numOfElements = ori -> numOfElements;
	ret -> elementList = malloc(sizeof(void *)*ret -> capacityOfElements);
	memcpy(ret -> elementList, ori -> elementList, sizeof(void *)*ret -> capacityOfElements);
	ret -> elemDeallocator = ori -> elemDeallocator;
	return ret;
}

void ArrayListDestroy(ArrayList * list){
	srInt_64 x1;
	if(list -> elemDeallocator)
		for(x1 = 0;x1 < list->numOfElements; x1++)
			list -> elemDeallocator(list -> elementList[x1]);
		
	free(list -> elementList);
	free(list);
}

void * ArrayListGet(ArrayList * list, srInt_64 n){
	if(n<0 || n >= list->numOfElements)return NULL;
	return list -> elementList[n];
}

int ArrayListPush_NoRepeatedPtr(ArrayList * list, void * new_elem){
	srInt_64 ii;
	for(ii = 0; ii < list->numOfElements; ii++){
		if(list -> elementList[ii] == new_elem) return -1;
	}
	return ArrayListPush(list,new_elem);
}

int ArrayListPush(ArrayList * list, void * new_elem){
	if(list -> capacityOfElements <= list->numOfElements){
		if(list -> capacityOfElements *1.3 > list -> capacityOfElements + 10)
				list -> capacityOfElements = list -> capacityOfElements *1.3;
		else	list -> capacityOfElements = list -> capacityOfElements + 10;

		list -> elementList=realloc(list -> elementList, sizeof(void *) * list -> capacityOfElements);
		assert(list -> elementList);
	}
	list->elementList[list->numOfElements++] = new_elem;
	return list->numOfElements;
}
void ArrayListSetDeallocationFunction(ArrayList * list,  void (*elem_deallocator)(void *elem)){
	list -> elemDeallocator = elem_deallocator;
}

int ArrayListSort_comp_pntr( void * L_elem, void * R_elem, ArrayList * me ){
	if( L_elem > R_elem )return 1;
	if( L_elem < R_elem )return -1;
	return 0;
}

int ArrayListSort_compare(void * sortdata0, int L, int R){
	void ** sortdata = sortdata0;
	ArrayList * list = sortdata[0];
	int (*comp_elems)(void * L_elem, void * R_elem, ArrayList * me) = sortdata[1];

	void * L_elem = list -> elementList[L];
	void * R_elem = list -> elementList[R];
	return comp_elems(L_elem, R_elem, list);
}

void ArrayListSort_exchange(void * sortdata0, int L, int R){
	void ** sortdata = sortdata0;
	ArrayList * list = sortdata[0];

	void * tmpp = list -> elementList[L];
	list -> elementList[L] = list -> elementList[R];
	list -> elementList[R] = tmpp;
}

void ArrayListSort_merge(void * sortdata0, int start, int items, int items2){
	void ** sortdata = sortdata0;
	ArrayList * list = sortdata[0];
	int (*comp_elems)(void * L_elem, void * R_elem, ArrayList * me) = sortdata[1];

	void ** merged = malloc(sizeof(void *)*(items + items2));
	int write_cursor, read1=start, read2=start+items;

	for(write_cursor = 0; write_cursor < items + items2; write_cursor++){
		void * Elm1 = (read1 == start + items)? NULL:list -> elementList[read1];
		void * Elm2 = (read2 == start + items + items2)?NULL:list -> elementList[read2];

		int select_1 = (read1 == start + items)?0:( read2 == start + items + items2 || comp_elems(Elm1, Elm2, list) < 0 );
		if(select_1) merged[write_cursor] = list -> elementList[read1++];
		else merged[write_cursor] = list -> elementList[read2++];
		if(read2 > start + items + items2){
			SUBREADprintf("R1: %d < %d ; R2: %d < %d\n", read1,  start + items, read2, start + items + items2);
		}
		assert(read2 <=start + items + items2);
		assert(read1 <=start + items);
	}

	memcpy(list -> elementList + start, merged, sizeof(void *) * (items + items2));
	free(merged);
}

void ArrayListSort(ArrayList * list, int compare_L_minus_R(void * L_elem, void * R_elem, ArrayList * me)){
	void * sortdata[2];	
	sortdata[0] = list;
	sortdata[1] = compare_L_minus_R?compare_L_minus_R:ArrayListSort_comp_pntr;

	if(list -> numOfElements > 20)
		merge_sort(sortdata, list -> numOfElements, ArrayListSort_compare, ArrayListSort_exchange, ArrayListSort_merge);
	else
		basic_sort(sortdata, list -> numOfElements, ArrayListSort_compare, ArrayListSort_exchange);
}

int ArrayListLLUComparison(void * L_elem, void * R_elem){
	srUInt_64 lint = L_elem-NULL;
	srUInt_64 rint = R_elem-NULL;
	if(lint<rint)return -1;
	if(lint>rint)return 1;
	return 0;
}

srInt_64 ArrayListFindNextDent(ArrayList * list, srUInt_64 value_less_than_dent){
	srInt_64 h=list->numOfElements-1,l=0,m=-1l;
	if( list -> elementList[h]- NULL <= value_less_than_dent )return -1l;
	while(h>l){
		m=(h+l)/2;
		srUInt_64 mv = list -> elementList[m] - NULL;
		if(mv < value_less_than_dent) l=m+1;
		else if(mv > value_less_than_dent) h=m-1;
		else break;
	}
	if(m<2)m=0;
	else m-=2;

	for(; m>=0 && list -> elementList[m] - NULL >= value_less_than_dent; m--);

	for(m = max(0, m); m < list->numOfElements ; m++){
		if( list -> elementList[m] - NULL > value_less_than_dent )return m;
	}
	SUBREADprintf("ALGORITHM ERROR!! DID YOU SORT THE LIST???\n");
	return -2l;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableCreate() - creates a new HashTable
 *  DESCRIPTION:
 *	  Creates a new HashTable.  When finished with this HashTable, it
 *	  should be explicitly destroyed by calling the HashTableDestroy()
 *	  function.
 *  EFFICIENCY:
 *	  O(1)
 *  ARGUMENTS:
 *	  numOfBuckets - the number of buckets to start the HashTable out with.
 *					 Must be greater than zero, and should be prime.
 *					 Ideally, the number of buckets should between 1/5
 *					 and 1 times the expected number of elements in the
 *					 HashTable.  Values much more or less than this will
 *					 result in wasted memory or decreased performance
 *					 respectively.  The number of buckets in a HashTable
 *					 can be re-calculated to an appropriate number by
 *					 calling the HashTableRehash() function once the
 *					 HashTable has been populated.  The number of buckets
 *					 in a HashTable may also be re-calculated
 *					 automatically if the ratio of elements to buckets
 *					 passes the thresholds set by HashTableSetIdealRatio().
 *  RETURNS:
 *	  HashTable	- a new Hashtable, or NULL on error
\*--------------------------------------------------------------------------*/

HashTable * StringTableCreate(srInt_64 numOfBuckets){
	HashTable * ret = HashTableCreate(numOfBuckets);
	HashTableSetHashFunction(ret, HashTableStringHashFunction);
	HashTableSetKeyComparisonFunction(ret , my_strcmp);
	return ret;
}

void StringTableReverse_Run(void * ky, void * va, HashTable * tab){
	HashTable * res = tab -> appendix1;
	HashTablePut(res, va, ky);
}

HashTable * StringTableReverse(HashTable * ori){
	HashTable * res = StringTableCreate(ori -> numOfBuckets);
	void * oriap1 = ori -> appendix1;
	ori -> appendix1 = res;

	HashTableIteration(ori, StringTableReverse_Run);
	ori -> appendix1 = oriap1;
	return res;
}

ArrayList * HashTableKeys(HashTable * tab){
	int i;
	ArrayList * ret = ArrayListCreate(tab -> numOfElements);
	for (i=0; i< tab ->numOfBuckets; i++) {
		KeyValuePair *pair = tab ->bucketArray[i];
		while (pair != NULL) {
			ArrayListPush(ret, pair -> key);
			KeyValuePair *nextPair = pair->next;
			pair = nextPair;
		}
	}
	return ret;
}

HashTable *HashTableCreate(srInt_64 numOfBuckets) {
	HashTable *hashTable;
	int i;


	assert(numOfBuckets > 0);

	hashTable = (HashTable *) malloc(sizeof(HashTable));
	if (hashTable == NULL)
		return NULL;

	hashTable->appendix1=NULL;
	hashTable->appendix2=NULL;
	hashTable->appendix3=NULL;

	hashTable->counter1=0;
	hashTable->counter2=0;
	hashTable->counter3=0;

	hashTable->bucketArray = (KeyValuePair **)
						malloc(numOfBuckets * sizeof(KeyValuePair *));
	if (hashTable->bucketArray == NULL) {
		free(hashTable);
		return NULL;
	}
	
	hashTable->numOfBuckets = numOfBuckets;
	hashTable->numOfElements = 0;

	for (i=0; i<numOfBuckets; i++)
		hashTable->bucketArray[i] = NULL;

	hashTable->idealRatio = 3.0;
	hashTable->lowerRehashThreshold = 0.0;
	hashTable->upperRehashThreshold = 15.0;

	hashTable->keycmp = pointercmp;
	hashTable->valuecmp = pointercmp;
	hashTable->hashFunction = pointerHashFunction;
	hashTable->keyDeallocator = NULL;
	hashTable->valueDeallocator = NULL;

	return hashTable;
}


/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableDestroy() - destroys an existing HashTable
 *  DESCRIPTION:
 *	  Destroys an existing HashTable.
 *  EFFICIENCY:
 *	  O(n)
 *  ARGUMENTS:
 *	  hashTable	- the HashTable to destroy
 *  RETURNS:
 *	  <nothing>
\*--------------------------------------------------------------------------*/

void HashTableDestroy(HashTable *hashTable) {
	int i;

	for (i=0; i<hashTable->numOfBuckets; i++) {
		KeyValuePair *pair = hashTable->bucketArray[i];
		while (pair != NULL) {
			KeyValuePair *nextPair = pair->next;

			if (hashTable->keyDeallocator != NULL)
				hashTable->keyDeallocator((void *) pair->key);
			if (hashTable->valueDeallocator != NULL){
//			fprintf(stderr,"FREE %p\n", pair->value);
			hashTable->valueDeallocator(pair->value);
		}
			free(pair);
			pair = nextPair;
		}
	}

	free(hashTable->bucketArray);
	free(hashTable);
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableContainsKey() - checks the existence of a key in a HashTable
 *  DESCRIPTION:
 *	  Determines whether or not the specified HashTable contains the
 *	  specified key.  Uses the comparison function specified by
 *	  HashTableSetKeyComparisonFunction().
 *  EFFICIENCY:
 *	  O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *	  hashTable	- the HashTable to search
 *	  key		  - the key to search for
 *  RETURNS:
 *	  bool		 - whether or not the specified HashTable contains the
 *					 specified key.
\*--------------------------------------------------------------------------*/

int HashTableContainsKey(const HashTable *hashTable, const void *key) {
	return (HashTableGet(hashTable, key) != NULL);
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableContainsValue()
 *						 - checks the existence of a value in a HashTable
 *  DESCRIPTION:
 *	  Determines whether or not the specified HashTable contains the
 *	  specified value.  Unlike HashTableContainsKey(), this function is
 *	  not very efficient, since it has to scan linearly looking for a
 *	  match.  Uses the comparison function specified by
 *	  HashTableSetValueComparisonFunction().
 *  EFFICIENCY:
 *	  O(n)
 *  ARGUMENTS:
 *	  hashTable	- the HashTable to search
 *	  value		- the value to search for
 *  RETURNS:
 *	  bool		 - whether or not the specified HashTable contains the
 *					 specified value.
\*--------------------------------------------------------------------------*/

int HashTableContainsValue(const HashTable *hashTable, const void *value) {
	int i;

	for (i=0; i<hashTable->numOfBuckets; i++) {
		KeyValuePair *pair = hashTable->bucketArray[i];
		while (pair != NULL) {
			if (hashTable->valuecmp(value, pair->value) == 0)
				return 1;
			pair = pair->next;
		}
	}

	return 0;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTablePut() - adds a key/value pair to a HashTable
 *  DESCRIPTION:
 *	  Adds the specified key/value pair to the specified HashTable.  If
 *	  the key already exists in the HashTable (determined by the comparison
 *	  function specified by HashTableSetKeyComparisonFunction()), its value
 *	  is replaced by the new value.  May trigger an auto-rehash (see
 *	  HashTableSetIdealRatio()).  It is illegal to specify NULL as the
 *	  key or value.
 *  EFFICIENCY:
 *	  O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *	  hashTable	- the HashTable to add to
 *	  key		  - the key to add or whose value to replace
 *	  value		- the value associated with the key
 *  RETURNS:
 *	  err		  - 0 if successful, -1 if an error was encountered
\*--------------------------------------------------------------------------*/

int HashTablePut(HashTable *hashTable, const void *key, void *value) {
	return HashTablePutReplace(hashTable, key, value, 1);
}

int HashTablePutReplaceEx(HashTable *hashTable, const void *key, void *value, int replace_key, int dealloc_key, int dealloc_value) {
	srInt_64 hashValue;
	KeyValuePair *pair;

	assert(key != NULL);
	assert(value != NULL);

	hashValue = hashTable->hashFunction(key) % hashTable->numOfBuckets;


	pair = hashTable->bucketArray[hashValue];

	while (pair != NULL && hashTable->keycmp(key, pair->key) != 0)
		pair = pair->next;

	if (pair) {
		if (pair->key != key) {
			if(replace_key) {
				if(hashTable->keyDeallocator && dealloc_key)
					hashTable->keyDeallocator((void *) pair->key);
				pair->key = (void *)key;
			}
		}
		if (pair->value != value) {
			if (hashTable->valueDeallocator != NULL && dealloc_value)
				hashTable->valueDeallocator(pair->value);
			pair->value = (void *)value;
		}
	}
	else {
		KeyValuePair *newPair = (KeyValuePair *) malloc(sizeof(KeyValuePair));
		if (newPair == NULL) {
		return -1;
		}
		else {
			newPair->key = (void *)key;
			newPair->value = (void *)value;
			newPair->next = hashTable->bucketArray[hashValue];
			hashTable->bucketArray[hashValue] = newPair;
			hashTable->numOfElements++;

			if (hashTable->upperRehashThreshold > hashTable->idealRatio) {
				float elementToBucketRatio = (float) hashTable->numOfElements /
											 (float) hashTable->numOfBuckets;
				if (elementToBucketRatio > hashTable->upperRehashThreshold)
					HashTableRehash(hashTable, 0);
			}
		}
	}

	return 0;
}
int HashTablePutReplace(HashTable *hashTable, const void *key, void *value, int replace_key) {
	return HashTablePutReplaceEx(hashTable, key, value,  replace_key, 1, 1);
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableGet() - retrieves the value of a key in a HashTable
 *  DESCRIPTION:
 *	  Retrieves the value of the specified key in the specified HashTable.
 *	  Uses the comparison function specified by
 *	  HashTableSetKeyComparisonFunction().
 *  EFFICIENCY:
 *	  O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *	  hashTable	- the HashTable to search
 *	  key		  - the key whose value is desired
 *  RETURNS:
 *	  void *	   - the value of the specified key, or NULL if the key
 *					 doesn't exist in the HashTable
\*--------------------------------------------------------------------------*/

void *HashTableGetKey(const HashTable *hashTable, const void *key) {
	srInt_64 hashValue = hashTable->hashFunction(key) % hashTable->numOfBuckets;

	KeyValuePair *pair = hashTable->bucketArray[hashValue];

	while (pair != NULL && hashTable->keycmp(key, pair->key) != 0)
		pair = pair->next;

	return (pair == NULL)? NULL : pair->key;
}


void *HashTableGet(const HashTable *hashTable, const void *key) {
	srInt_64 hashValue = hashTable->hashFunction(key) % hashTable->numOfBuckets;

	KeyValuePair *pair = hashTable->bucketArray[hashValue];

	while (pair != NULL && hashTable->keycmp(key, pair->key) != 0)
		pair = pair->next;

	return (pair == NULL)? NULL : pair->value;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableRemove() - removes a key/value pair from a HashTable
 *  DESCRIPTION:
 *	  Removes the key/value pair identified by the specified key from the
 *	  specified HashTable if the key exists in the HashTable.  May trigger
 *	  an auto-rehash (see HashTableSetIdealRatio()).
 *  EFFICIENCY:
 *	  O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *	  hashTable	- the HashTable to remove the key/value pair from
 *	  key		  - the key specifying the key/value pair to be removed
 *  RETURNS:
 *	  <nothing>
\*--------------------------------------------------------------------------*/

void HashTableRemove(HashTable *hashTable, const void *key) {
	srInt_64 hashValue = hashTable->hashFunction(key) % hashTable->numOfBuckets;


	KeyValuePair *pair = hashTable->bucketArray[hashValue];
	KeyValuePair *previousPair = NULL;

	while (pair != NULL && hashTable->keycmp(key, pair->key) != 0) {
		previousPair = pair;
		pair = pair->next;
	}

	if (pair != NULL) {
		if (hashTable->keyDeallocator != NULL)
			hashTable->keyDeallocator((void *) pair->key);
		if (hashTable->valueDeallocator != NULL)
			hashTable->valueDeallocator(pair->value);
		if (previousPair != NULL)
			previousPair->next = pair->next;
		else
			hashTable->bucketArray[hashValue] = pair->next;
		free(pair);
		hashTable->numOfElements--;

		if (hashTable->lowerRehashThreshold > 0.0) {
			float elementToBucketRatio = (float) hashTable->numOfElements /
										 (float) hashTable->numOfBuckets;
			if (elementToBucketRatio < hashTable->lowerRehashThreshold)
				HashTableRehash(hashTable, 0);
		}
	}


}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableRemoveAll() - removes all key/value pairs from a HashTable
 *  DESCRIPTION:
 *	  Removes all key/value pairs from the specified HashTable.  May trigger
 *	  an auto-rehash (see HashTableSetIdealRatio()).
 *  EFFICIENCY:
 *	  O(n)
 *  ARGUMENTS:
 *	  hashTable	- the HashTable to remove all key/value pairs from
 *  RETURNS:
 *	  <nothing>
\*--------------------------------------------------------------------------*/

void HashTableRemoveAll(HashTable *hashTable) {
	int i;

	for (i=0; i<hashTable->numOfBuckets; i++) {
		KeyValuePair *pair = hashTable->bucketArray[i];
		while (pair != NULL) {
			KeyValuePair *nextPair = pair->next;
			if (hashTable->keyDeallocator != NULL)
				hashTable->keyDeallocator((void *) pair->key);
			if (hashTable->valueDeallocator != NULL)
				hashTable->valueDeallocator(pair->value);
			free(pair);
			pair = nextPair;
		}
		hashTable->bucketArray[i] = NULL;
	}

	hashTable->numOfElements = 0;
	HashTableRehash(hashTable, 5);
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableIsEmpty() - determines if a HashTable is empty
 *  DESCRIPTION:
 *	  Determines whether or not the specified HashTable contains any
 *	  key/value pairs.
 *  EFFICIENCY:
 *	  O(1)
 *  ARGUMENTS:
 *	  hashTable	- the HashTable to check
 *  RETURNS:
 *	  bool		 - whether or not the specified HashTable contains any
 *					 key/value pairs
\*--------------------------------------------------------------------------*/

int HashTableIsEmpty(const HashTable *hashTable) {
	return (hashTable->numOfElements == 0);
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableSize() - returns the number of elements in a HashTable
 *  DESCRIPTION:
 *	  Returns the number of key/value pairs that are present in the
 *	  specified HashTable.
 *  EFFICIENCY:
 *	  O(1)
 *  ARGUMENTS:
 *	  hashTable	- the HashTable whose size is requested
 *  RETURNS:
 *	  long		 - the number of key/value pairs that are present in
 *					 the specified HashTable
\*--------------------------------------------------------------------------*/

srInt_64 HashTableSize(const HashTable *hashTable) {
	return hashTable->numOfElements;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableGetNumBuckets() - returns the number of buckets in a HashTable
 *  DESCRIPTION:
 *	  Returns the number of buckets that are in the specified HashTable.
 *	  This may change dynamically throughout the life of a HashTable if
 *	  automatic or manual rehashing is performed.
 *  EFFICIENCY:
 *	  O(1)
 *  ARGUMENTS:
 *	  hashTable	- the HashTable whose number of buckets is requested
 *  RETURNS:
 *	  long		 - the number of buckets that are in the specified
 *					 HashTable
\*--------------------------------------------------------------------------*/

srInt_64 HashTableGetNumBuckets(const HashTable *hashTable) {
	return hashTable->numOfBuckets;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableSetKeyComparisonFunction()
 *			  - specifies the function used to compare keys in a HashTable
 *  DESCRIPTION:
 *	  Specifies the function used to compare keys in the specified
 *	  HashTable.  The specified function should return zero if the two
 *	  keys are considered equal, and non-zero otherwise.  The default
 *	  function is one that simply compares pointers.
 *  ARGUMENTS:
 *	  hashTable	- the HashTable whose key comparison function is being
 *					 specified
 *	  keycmp	   - a function which returns zero if the two arguments
 *					 passed to it are considered "equal" keys and non-zero
 *					 otherwise
 *  RETURNS:
 *	  <nothing>
\*--------------------------------------------------------------------------*/

void HashTableSetKeyComparisonFunction(HashTable *hashTable,
		int (*keycmp)(const void *key1, const void *key2)) {
	assert(keycmp != NULL);
	hashTable->keycmp = keycmp;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableSetValueComparisonFunction()
 *			  - specifies the function used to compare values in a HashTable
 *  DESCRIPTION:
 *	  Specifies the function used to compare values in the specified
 *	  HashTable.  The specified function should return zero if the two
 *	  values are considered equal, and non-zero otherwise.  The default
 *	  function is one that simply compares pointers.
 *  ARGUMENTS:
 *	  hashTable	- the HashTable whose value comparison function is being
 *					 specified
 *	  valuecmp	 - a function which returns zero if the two arguments
 *					 passed to it are considered "equal" values and non-zero
 *					 otherwise
 *  RETURNS:
 *	  <nothing>
\*--------------------------------------------------------------------------*/

void HashTableSetValueComparisonFunction(HashTable *hashTable,
		int (*valuecmp)(const void *value1, const void *value2)) {
	assert(valuecmp != NULL);
	hashTable->valuecmp = valuecmp;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableSetHashFunction()
 *			  - specifies the hash function used by a HashTable
 *  DESCRIPTION:
 *	  Specifies the function used to determine the hash value for a key
 *	  in the specified HashTable (before modulation).  An ideal hash
 *	  function is one which is easy to compute and approximates a
 *	  "random" function.  The default function is one that works
 *	  relatively well for pointers.  If the HashTable keys are to be
 *	  strings (which is probably the case), then this default function
 *	  will not suffice, in which case consider using the provided
 *	  HashTableStringHashFunction() function.
 *  ARGUMENTS:
 *	  hashTable	- the HashTable whose hash function is being specified
 *	  hashFunction - a function which returns an appropriate hash code
 *					 for a given key
 *  RETURNS:
 *	  <nothing>
\*--------------------------------------------------------------------------*/

void HashTableSetHashFunction(HashTable *hashTable,
		srUInt_64 (*hashFunction)(const void *key))
{
	assert(hashFunction != NULL);
	hashTable->hashFunction = hashFunction;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableRehash() - reorganizes a HashTable to be more efficient
 *  DESCRIPTION:
 *	  Reorganizes a HashTable to be more efficient.  If a number of
 *	  buckets is specified, the HashTable is rehashed to that number of
 *	  buckets.  If 0 is specified, the HashTable is rehashed to a number
 *	  of buckets which is automatically calculated to be a prime number
 *	  that achieves (as closely as possible) the ideal element-to-bucket 
 *	  ratio specified by the HashTableSetIdealRatio() function.
 *  EFFICIENCY:
 *	  O(n)
 *  ARGUMENTS:
 *	  hashTable	- the HashTable to be reorganized
 *	  numOfBuckets - the number of buckets to rehash the HashTable to.
 *					 Should be prime.  Ideally, the number of buckets
 *					 should be between 1/5 and 1 times the expected
 *					 number of elements in the HashTable.  Values much
 *					 more or less than this will result in wasted memory
 *					 or decreased performance respectively.  If 0 is
 *					 specified, an appropriate number of buckets is
 *					 automatically calculated.
 *  RETURNS:
 *	  <nothing>
\*--------------------------------------------------------------------------*/

void HashTableRehash(HashTable *hashTable, srInt_64 numOfBuckets) {
	KeyValuePair **newBucketArray;
	int i;

	assert(numOfBuckets >= 0);
	if (numOfBuckets == 0)
		numOfBuckets = calculateIdealNumOfBuckets(hashTable);

	if (numOfBuckets == hashTable->numOfBuckets)
		return; /* already the right size! */

	newBucketArray = (KeyValuePair **)
								malloc(numOfBuckets * sizeof(KeyValuePair *));
	if (newBucketArray == NULL) {
		/* Couldn't allocate memory for the new array.  This isn't a fatal
		 * error; we just can't perform the rehash. */
		return;
	}

	for (i=0; i<numOfBuckets; i++)
		newBucketArray[i] = NULL;

	for (i=0; i<hashTable->numOfBuckets; i++) {
		KeyValuePair *pair = hashTable->bucketArray[i];
		while (pair != NULL) {
			KeyValuePair *nextPair = pair->next;
			srInt_64 hashValue = hashTable->hashFunction(pair->key) % numOfBuckets;
			pair->next = newBucketArray[hashValue];
			newBucketArray[hashValue] = pair;
			pair = nextPair;
		}
	}

	free(hashTable->bucketArray);
	hashTable->bucketArray = newBucketArray;
	hashTable->numOfBuckets = numOfBuckets;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableSetIdealRatio()
 *			  - sets the ideal element-to-bucket ratio of a HashTable
 *  DESCRIPTION:
 *	  Sets the ideal element-to-bucket ratio, as well as the lower and
 *	  upper auto-rehash thresholds, of the specified HashTable.  Note
 *	  that this function doesn't actually perform a rehash.
 *
 *	  The default values for these properties are 3.0, 0.0 and 15.0
 *	  respectively.  This is likely fine for most situations, so there
 *	  is probably no need to call this function.
 *  ARGUMENTS:
 *	  hashTable	- a HashTable
 *	  idealRatio   - the ideal element-to-bucket ratio.  When a rehash
 *					 occurs (either manually via a call to the
 *					 HashTableRehash() function or automatically due the
 *					 the triggering of one of the thresholds below), the
 *					 number of buckets in the HashTable will be
 *					 recalculated to be a prime number that achieves (as
 *					 closely as possible) this ideal ratio.  Must be a
 *					 positive number.
 *	  lowerRehashThreshold
 *				   - the element-to-bucket ratio that is considered
 *					 unacceptably low (i.e., too few elements per bucket).
 *					 If the actual ratio falls below this number, a
 *					 rehash will automatically be performed.  Must be
 *					 lower than the value of idealRatio.  If no ratio
 *					 is considered unacceptably low, a value of 0.0 can
 *					 be specified.
 *	  upperRehashThreshold
 *				   - the element-to-bucket ratio that is considered
 *					 unacceptably high (i.e., too many elements per bucket).
 *					 If the actual ratio rises above this number, a
 *					 rehash will automatically be performed.  Must be
 *					 higher than idealRatio.  However, if no ratio
 *					 is considered unacceptably high, a value of 0.0 can
 *					 be specified.
 *  RETURNS:
 *	  <nothing>
\*--------------------------------------------------------------------------*/

void HashTableSetIdealRatio(HashTable *hashTable, float idealRatio,
		float lowerRehashThreshold, float upperRehashThreshold) {
	assert(idealRatio > 0.0);
	assert(lowerRehashThreshold < idealRatio);
	assert(upperRehashThreshold == 0.0 || upperRehashThreshold > idealRatio);

	hashTable->idealRatio = idealRatio;
	hashTable->lowerRehashThreshold = lowerRehashThreshold;
	hashTable->upperRehashThreshold = upperRehashThreshold;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableSetDeallocationFunctions()
 *			  - sets the key and value deallocation functions of a HashTable
 *  DESCRIPTION:
 *	  Sets the key and value deallocation functions of the specified
 *	  HashTable.  This determines what happens to a key or a value when it
 *	  is removed from the HashTable.  If the deallocation function is NULL
 *	  (the default if this function is never called), its reference is
 *	  simply dropped and it is up to the calling program to perform the
 *	  proper memory management.  If the deallocation function is non-NULL,
 *	  it is called to free the memory used by the object.  E.g., for simple
 *	  objects, an appropriate deallocation function may be free().
 *
 *	  This affects the behaviour of the HashTableDestroy(), HashTablePut(),
 *	  HashTableRemove() and HashTableRemoveAll() functions.
 *  ARGUMENTS:
 *	  hashTable	- a HashTable
 *	  keyDeallocator
 *				   - if non-NULL, the function to be called when a key is
 *					 removed from the HashTable.
 *	  valueDeallocator
 *				   - if non-NULL, the function to be called when a value is
 *					 removed from the HashTable.
 *  RETURNS:
 *	  <nothing>
\*--------------------------------------------------------------------------*/

void HashTableSetDeallocationFunctions(HashTable *hashTable,
		void (*keyDeallocator)(void *key),
		void (*valueDeallocator)(void *value)) {
	hashTable->keyDeallocator = keyDeallocator;
	hashTable->valueDeallocator = valueDeallocator;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *	  HashTableStringHashFunction() - a good hash function for strings
 *  DESCRIPTION:
 *	  A hash function that is appropriate for hashing strings.  Note that
 *	  this is not the default hash function.  To make it the default hash
 *	  function, call HashTableSetHashFunction(hashTable,
 *	  HashTableStringHashFunction).
 *  ARGUMENTS:
 *	  key	- the key to be hashed
 *  RETURNS:
 *	  unsigned long - the unmodulated hash value of the key
\*--------------------------------------------------------------------------*/

srUInt_64 HashTableStringHashFunction(const void *key) {
	const unsigned char *str = (const unsigned char *) key;
	srUInt_64 hashValue = 0;
	int i;

	for (i=0; str[i] != '\0'; i++)
		hashValue = hashValue * 37 + str[i];

	return hashValue;
}

static int pointercmp(const void *pointer1, const void *pointer2) {
	return (pointer1 != pointer2);
}

static srUInt_64 pointerHashFunction(const void *pointer) {
	return (srUInt_64)(pointer - NULL);
}

static int isProbablePrime(srInt_64 oddNumber) {
	srInt_64 i;

	for (i=3; i<51; i+=2)
		if (oddNumber == i)
			return 1;
		else if (oddNumber%i == 0)
			return 0;

	return 1; /* maybe */
}

static srInt_64 calculateIdealNumOfBuckets(HashTable *hashTable) {
	srInt_64 idealNumOfBuckets = hashTable->numOfElements / hashTable->idealRatio;
	if (idealNumOfBuckets < 5)
		idealNumOfBuckets = 5;
	else
		idealNumOfBuckets |= 0x01; /* make it an odd number */
	while (!isProbablePrime(idealNumOfBuckets))
		idealNumOfBuckets += 2;

	return idealNumOfBuckets;
}


void free_values_destroy(HashTable * tab)
{

	KeyValuePair * cursor;
	int bucket;
	
	for(bucket=0; bucket< tab -> numOfBuckets; bucket++)
	{
		cursor = tab -> bucketArray[bucket];
		while (1)
		{
			if(!cursor) break;
			char * read_txt = (char *) cursor ->value;
			free(read_txt);
			cursor = cursor->next;
		}
	}

	HashTableDestroy(tab);
}

ArrayList * HashTableKeyArray(HashTable * tab){
	ArrayList *ret = ArrayListCreate(20);
	int i;
	for (i=0; i< tab ->numOfBuckets; i++) {
		KeyValuePair *pair = tab ->bucketArray[i];
		while (pair != NULL) {
			ArrayListPush(ret, (void *)pair -> key);
			KeyValuePair *nextPair = pair->next;
			pair = nextPair;
		}
	}
	return ret;
}

void HashTableIteration(HashTable * tab, void process_item(void * key, void * hashed_obj, HashTable * tab) )
{
	int i;
	for (i=0; i< tab ->numOfBuckets; i++) {
		KeyValuePair *pair = tab ->bucketArray[i];
		while (pair != NULL) {
			process_item(( void * )pair -> key, pair -> value, tab);
			KeyValuePair *nextPair = pair->next;
			pair = nextPair;
		}
	}
}

void HashTableSortedIndexes_copy_idx( void *k, void *v, HashTable * k2int_tab ){
	ArrayList * dsta = k2int_tab -> appendix1;
	ArrayListPush(dsta, k);
}

int HashTableSortedIndexes_cmp_idx( void * Lv, void * Rv, ArrayList * me ){
	void ** appdx = me -> appendix1;
	HashTable * k2int_tab = appdx[0];
	void * Lhash = HashTableGet(k2int_tab, Lv);
	void * Rhash = HashTableGet(k2int_tab, Rv);
	void * large_first = appdx[1];
	if(large_first){
		if(Lhash > Rhash) return -1;
		if(Lhash < Rhash) return 1;
	}else{
		if(Lhash > Rhash) return 1;
		if(Lhash < Rhash) return -1;
	}
	return 0;
}

ArrayList * HashTableSortedIndexes(HashTable * k2int_tab, int larger_value_first){
	ArrayList *ret = HashTableKeys(k2int_tab);
	void * appx[2];
	ret -> appendix1 = appx;
	appx[0]=k2int_tab;
	appx[1]=NULL+larger_value_first;
	ArrayListSort(ret, HashTableSortedIndexes_cmp_idx);
	return ret;
}

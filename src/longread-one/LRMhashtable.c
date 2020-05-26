/*--------------------------------------------------------------------------*\
 *                   -----===== HashTable =====-----
 *
 * Author: Keith Pomakis (pomakis@pobox.com)
 * Date:   August, 1998
 * Released to the public domain.
 *
 *--------------------------------------------------------------------------
 * $Id: LRMhashtable.c,v 1.5 2019/08/28 22:31:15 cvs Exp $
\*--------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include "LRMconfig.h"
#include "LRMhashtable.h"

static int LRMpointercmp(const void *pointer1, const void *pointer2);
static srUInt_64 LRMpointerHashFunction(const void *pointer);
static int LRMisProbablePrime(srInt_64 number);
static srInt_64 LRMcalculateIdealNumOfBuckets(HashTable *hashTable);


ArrayList * LRMArrayListCreate(int init_capacity){
	ArrayList * ret = malloc(sizeof(ArrayList));
	memset(ret,0,sizeof(ArrayList));
	ret -> capacityOfElements = init_capacity;
	ret -> elementList = malloc(sizeof(void *)*init_capacity);
	return ret;
}

void LRMArrayListDestroy(ArrayList * list){
	srInt_64 x1;
	if(list -> elemDeallocator)
		for(x1 = 0;x1 < list->numOfElements; x1++)
			list -> elemDeallocator(list -> elementList[x1]);
		
	free(list -> elementList);
	free(list);
}

void * LRMArrayListGet(ArrayList * list, srInt_64 n){
	if(n<0 || n >= list->numOfElements)return NULL;
	return list -> elementList[n];
}

int LRMArrayListPush(ArrayList * list, void * new_elem){
	if(list -> capacityOfElements <= list->numOfElements){
		list -> capacityOfElements *=1.3;
		list -> elementList=realloc(list -> elementList, sizeof(void *)*list -> capacityOfElements);
	}
	list->elementList[list->numOfElements++] = new_elem;
	return list->numOfElements;
}
void LRMArrayListSetDeallocationFunction(ArrayList * list,  void (*elem_deallocator)(void *elem)){
	list -> elemDeallocator = elem_deallocator;
}



/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableCreate() - creates a new HashTable
 *  DESCRIPTION:
 *      Creates a new HashTable.  When finished with this HashTable, it
 *      should be explicitly destroyed by calling the HashTableDestroy()
 *      function.
 *  EFFICIENCY:
 *      O(1)
 *  ARGUMENTS:
 *      numOfBuckets - the number of buckets to start the HashTable out with.
 *                     Must be greater than zero, and should be prime.
 *                     Ideally, the number of buckets should between 1/5
 *                     and 1 times the expected number of elements in the
 *                     HashTable.  Values much more or less than this will
 *                     result in wasted memory or decreased performance
 *                     respectively.  The number of buckets in a HashTable
 *                     can be re-calculated to an appropriate number by
 *                     calling the HashTableRehash() function once the
 *                     HashTable has been populated.  The number of buckets
 *                     in a HashTable may also be re-calculated
 *                     automatically if the ratio of elements to buckets
 *                     passes the thresholds set by HashTableSetIdealRatio().
 *  RETURNS:
 *      HashTable    - a new Hashtable, or NULL on error
\*--------------------------------------------------------------------------*/


HashTable *LRMHashTableCreate(srInt_64 numOfBuckets) {
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

    hashTable->keycmp = LRMpointercmp;
    hashTable->valuecmp = LRMpointercmp;
    hashTable->hashFunction = LRMpointerHashFunction;
    hashTable->keyDeallocator = NULL;
    hashTable->valueDeallocator = NULL;

    return hashTable;
}


/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableDestroy() - destroys an existing HashTable
 *  DESCRIPTION:
 *      Destroys an existing HashTable.
 *  EFFICIENCY:
 *      O(n)
 *  ARGUMENTS:
 *      hashTable    - the HashTable to destroy
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableDestroy(HashTable *hashTable) {
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
 *      HashTableContainsKey() - checks the existence of a key in a HashTable
 *  DESCRIPTION:
 *      Determines whether or not the specified HashTable contains the
 *      specified key.  Uses the comparison function specified by
 *      HashTableSetKeyComparisonFunction().
 *  EFFICIENCY:
 *      O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *      hashTable    - the HashTable to search
 *      key          - the key to search for
 *  RETURNS:
 *      bool         - whether or not the specified HashTable contains the
 *                     specified key.
\*--------------------------------------------------------------------------*/

int LRMHashTableContainsKey(const HashTable *hashTable, const void *key) {
    return (LRMHashTableGet(hashTable, key) != NULL);
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableContainsValue()
 *                         - checks the existence of a value in a HashTable
 *  DESCRIPTION:
 *      Determines whether or not the specified HashTable contains the
 *      specified value.  Unlike HashTableContainsKey(), this function is
 *      not very efficient, since it has to scan linearly looking for a
 *      match.  Uses the comparison function specified by
 *      HashTableSetValueComparisonFunction().
 *  EFFICIENCY:
 *      O(n)
 *  ARGUMENTS:
 *      hashTable    - the HashTable to search
 *      value        - the value to search for
 *  RETURNS:
 *      bool         - whether or not the specified HashTable contains the
 *                     specified value.
\*--------------------------------------------------------------------------*/

int LRMHashTableContainsValue(const HashTable *hashTable, const void *value) {
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
 *      HashTablePut() - adds a key/value pair to a HashTable
 *  DESCRIPTION:
 *      Adds the specified key/value pair to the specified HashTable.  If
 *      the key already exists in the HashTable (determined by the comparison
 *      function specified by HashTableSetKeyComparisonFunction()), its value
 *      is replaced by the new value.  May trigger an auto-rehash (see
 *      HashTableSetIdealRatio()).  It is illegal to specify NULL as the
 *      key or value.
 *  EFFICIENCY:
 *      O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *      hashTable    - the HashTable to add to
 *      key          - the key to add or whose value to replace
 *      value        - the value associated with the key
 *  RETURNS:
 *      err          - 0 if successful, -1 if an error was encountered
\*--------------------------------------------------------------------------*/

int LRMHashTablePut(HashTable *hashTable, const void *key, void *value) {
	return LRMHashTablePutReplace(hashTable, key, value, 1);
}
int LRMHashTablePutReplace(HashTable *hashTable, const void *key, void *value, int replace_key) {
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
            if(replace_key)
	    {
		if(hashTable->keyDeallocator)
			hashTable->keyDeallocator((void *) pair->key);
                pair->key = key;
	    }
        }
        if (pair->value != value) {
            if (hashTable->valueDeallocator != NULL)
                hashTable->valueDeallocator(pair->value);
            pair->value = value;
        }
    }
    else {
        KeyValuePair *newPair = (KeyValuePair *) malloc(sizeof(KeyValuePair));
        if (newPair == NULL) {
		return -1;
        }
        else {
            newPair->key = key;
            newPair->value = value;
            newPair->next = hashTable->bucketArray[hashValue];
            hashTable->bucketArray[hashValue] = newPair;
            hashTable->numOfElements++;

            if (hashTable->upperRehashThreshold > hashTable->idealRatio) {
                float elementToBucketRatio = (float) hashTable->numOfElements /
                                             (float) hashTable->numOfBuckets;
                if (elementToBucketRatio > hashTable->upperRehashThreshold)
                    LRMHashTableRehash(hashTable, 0);
            }
        }
    }

    return 0;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableGet() - retrieves the value of a key in a HashTable
 *  DESCRIPTION:
 *      Retrieves the value of the specified key in the specified HashTable.
 *      Uses the comparison function specified by
 *      HashTableSetKeyComparisonFunction().
 *  EFFICIENCY:
 *      O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *      hashTable    - the HashTable to search
 *      key          - the key whose value is desired
 *  RETURNS:
 *      void *       - the value of the specified key, or NULL if the key
 *                     doesn't exist in the HashTable
\*--------------------------------------------------------------------------*/

void *LRMHashTableGet(const HashTable *hashTable, const void *key) {
    srInt_64 hashValue = hashTable->hashFunction(key) % hashTable->numOfBuckets;

    KeyValuePair *pair = hashTable->bucketArray[hashValue];

    while (pair != NULL && hashTable->keycmp(key, pair->key) != 0)
        pair = pair->next;

    return (pair == NULL)? NULL : pair->value;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableRemove() - removes a key/value pair from a HashTable
 *  DESCRIPTION:
 *      Removes the key/value pair identified by the specified key from the
 *      specified HashTable if the key exists in the HashTable.  May trigger
 *      an auto-rehash (see HashTableSetIdealRatio()).
 *  EFFICIENCY:
 *      O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *      hashTable    - the HashTable to remove the key/value pair from
 *      key          - the key specifying the key/value pair to be removed
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableRemove(HashTable *hashTable, const void *key) {
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
                LRMHashTableRehash(hashTable, 0);
        }
    }


}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableRemoveAll() - removes all key/value pairs from a HashTable
 *  DESCRIPTION:
 *      Removes all key/value pairs from the specified HashTable.  May trigger
 *      an auto-rehash (see HashTableSetIdealRatio()).
 *  EFFICIENCY:
 *      O(n)
 *  ARGUMENTS:
 *      hashTable    - the HashTable to remove all key/value pairs from
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableRemoveAll(HashTable *hashTable) {
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
    LRMHashTableRehash(hashTable, 5);
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableIsEmpty() - determines if a HashTable is empty
 *  DESCRIPTION:
 *      Determines whether or not the specified HashTable contains any
 *      key/value pairs.
 *  EFFICIENCY:
 *      O(1)
 *  ARGUMENTS:
 *      hashTable    - the HashTable to check
 *  RETURNS:
 *      bool         - whether or not the specified HashTable contains any
 *                     key/value pairs
\*--------------------------------------------------------------------------*/

int LRMHashTableIsEmpty(const HashTable *hashTable) {
    return (hashTable->numOfElements == 0);
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableSize() - returns the number of elements in a HashTable
 *  DESCRIPTION:
 *      Returns the number of key/value pairs that are present in the
 *      specified HashTable.
 *  EFFICIENCY:
 *      O(1)
 *  ARGUMENTS:
 *      hashTable    - the HashTable whose size is requested
 *  RETURNS:
 *      long         - the number of key/value pairs that are present in
 *                     the specified HashTable
\*--------------------------------------------------------------------------*/

srInt_64 LRMHashTableSize(const HashTable *hashTable) {
    return hashTable->numOfElements;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableGetNumBuckets() - returns the number of buckets in a HashTable
 *  DESCRIPTION:
 *      Returns the number of buckets that are in the specified HashTable.
 *      This may change dynamically throughout the life of a HashTable if
 *      automatic or manual rehashing is performed.
 *  EFFICIENCY:
 *      O(1)
 *  ARGUMENTS:
 *      hashTable    - the HashTable whose number of buckets is requested
 *  RETURNS:
 *      long         - the number of buckets that are in the specified
 *                     HashTable
\*--------------------------------------------------------------------------*/

srInt_64 LRMHashTableGetNumBuckets(const HashTable *hashTable) {
    return hashTable->numOfBuckets;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableSetKeyComparisonFunction()
 *              - specifies the function used to compare keys in a HashTable
 *  DESCRIPTION:
 *      Specifies the function used to compare keys in the specified
 *      HashTable.  The specified function should return zero if the two
 *      keys are considered equal, and non-zero otherwise.  The default
 *      function is one that simply compares pointers.
 *  ARGUMENTS:
 *      hashTable    - the HashTable whose key comparison function is being
 *                     specified
 *      keycmp       - a function which returns zero if the two arguments
 *                     passed to it are considered "equal" keys and non-zero
 *                     otherwise
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableSetKeyComparisonFunction(HashTable *hashTable,
        int (*keycmp)(const void *key1, const void *key2)) {
    assert(keycmp != NULL);
    hashTable->keycmp = keycmp;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableSetValueComparisonFunction()
 *              - specifies the function used to compare values in a HashTable
 *  DESCRIPTION:
 *      Specifies the function used to compare values in the specified
 *      HashTable.  The specified function should return zero if the two
 *      values are considered equal, and non-zero otherwise.  The default
 *      function is one that simply compares pointers.
 *  ARGUMENTS:
 *      hashTable    - the HashTable whose value comparison function is being
 *                     specified
 *      valuecmp     - a function which returns zero if the two arguments
 *                     passed to it are considered "equal" values and non-zero
 *                     otherwise
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableSetValueComparisonFunction(HashTable *hashTable,
        int (*valuecmp)(const void *value1, const void *value2)) {
    assert(valuecmp != NULL);
    hashTable->valuecmp = valuecmp;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableSetHashFunction()
 *              - specifies the hash function used by a HashTable
 *  DESCRIPTION:
 *      Specifies the function used to determine the hash value for a key
 *      in the specified HashTable (before modulation).  An ideal hash
 *      function is one which is easy to compute and approximates a
 *      "random" function.  The default function is one that works
 *      relatively well for pointers.  If the HashTable keys are to be
 *      strings (which is probably the case), then this default function
 *      will not suffice, in which case consider using the provided
 *      HashTableStringHashFunction() function.
 *  ARGUMENTS:
 *      hashTable    - the HashTable whose hash function is being specified
 *      hashFunction - a function which returns an appropriate hash code
 *                     for a given key
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableSetHashFunction(HashTable *hashTable,
        srUInt_64 (*hashFunction)(const void *key))
{
    assert(hashFunction != NULL);
    hashTable->hashFunction = hashFunction;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableRehash() - reorganizes a HashTable to be more efficient
 *  DESCRIPTION:
 *      Reorganizes a HashTable to be more efficient.  If a number of
 *      buckets is specified, the HashTable is rehashed to that number of
 *      buckets.  If 0 is specified, the HashTable is rehashed to a number
 *      of buckets which is automatically calculated to be a prime number
 *      that achieves (as closely as possible) the ideal element-to-bucket 
 *      ratio specified by the HashTableSetIdealRatio() function.
 *  EFFICIENCY:
 *      O(n)
 *  ARGUMENTS:
 *      hashTable    - the HashTable to be reorganized
 *      numOfBuckets - the number of buckets to rehash the HashTable to.
 *                     Should be prime.  Ideally, the number of buckets
 *                     should be between 1/5 and 1 times the expected
 *                     number of elements in the HashTable.  Values much
 *                     more or less than this will result in wasted memory
 *                     or decreased performance respectively.  If 0 is
 *                     specified, an appropriate number of buckets is
 *                     automatically calculated.
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableRehash(HashTable *hashTable, srInt_64 numOfBuckets) {
    KeyValuePair **newBucketArray;
    int i;

    assert(numOfBuckets >= 0);
    if (numOfBuckets == 0)
        numOfBuckets = LRMcalculateIdealNumOfBuckets(hashTable);

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
 *      HashTableSetIdealRatio()
 *              - sets the ideal element-to-bucket ratio of a HashTable
 *  DESCRIPTION:
 *      Sets the ideal element-to-bucket ratio, as well as the lower and
 *      upper auto-rehash thresholds, of the specified HashTable.  Note
 *      that this function doesn't actually perform a rehash.
 *
 *      The default values for these properties are 3.0, 0.0 and 15.0
 *      respectively.  This is likely fine for most situations, so there
 *      is probably no need to call this function.
 *  ARGUMENTS:
 *      hashTable    - a HashTable
 *      idealRatio   - the ideal element-to-bucket ratio.  When a rehash
 *                     occurs (either manually via a call to the
 *                     HashTableRehash() function or automatically due the
 *                     the triggering of one of the thresholds below), the
 *                     number of buckets in the HashTable will be
 *                     recalculated to be a prime number that achieves (as
 *                     closely as possible) this ideal ratio.  Must be a
 *                     positive number.
 *      lowerRehashThreshold
 *                   - the element-to-bucket ratio that is considered
 *                     unacceptably low (i.e., too few elements per bucket).
 *                     If the actual ratio falls below this number, a
 *                     rehash will automatically be performed.  Must be
 *                     lower than the value of idealRatio.  If no ratio
 *                     is considered unacceptably low, a value of 0.0 can
 *                     be specified.
 *      upperRehashThreshold
 *                   - the element-to-bucket ratio that is considered
 *                     unacceptably high (i.e., too many elements per bucket).
 *                     If the actual ratio rises above this number, a
 *                     rehash will automatically be performed.  Must be
 *                     higher than idealRatio.  However, if no ratio
 *                     is considered unacceptably high, a value of 0.0 can
 *                     be specified.
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableSetIdealRatio(HashTable *hashTable, float idealRatio,
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
 *      HashTableSetDeallocationFunctions()
 *              - sets the key and value deallocation functions of a HashTable
 *  DESCRIPTION:
 *      Sets the key and value deallocation functions of the specified
 *      HashTable.  This determines what happens to a key or a value when it
 *      is removed from the HashTable.  If the deallocation function is NULL
 *      (the default if this function is never called), its reference is
 *      simply dropped and it is up to the calling program to perform the
 *      proper memory management.  If the deallocation function is non-NULL,
 *      it is called to free the memory used by the object.  E.g., for simple
 *      objects, an appropriate deallocation function may be free().
 *
 *      This affects the behaviour of the HashTableDestroy(), HashTablePut(),
 *      HashTableRemove() and HashTableRemoveAll() functions.
 *  ARGUMENTS:
 *      hashTable    - a HashTable
 *      keyDeallocator
 *                   - if non-NULL, the function to be called when a key is
 *                     removed from the HashTable.
 *      valueDeallocator
 *                   - if non-NULL, the function to be called when a value is
 *                     removed from the HashTable.
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableSetDeallocationFunctions(HashTable *hashTable,
        void (*keyDeallocator)(void *key),
        void (*valueDeallocator)(void *value)) {
    hashTable->keyDeallocator = keyDeallocator;
    hashTable->valueDeallocator = valueDeallocator;
}

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      HashTableStringHashFunction() - a good hash function for strings
 *  DESCRIPTION:
 *      A hash function that is appropriate for hashing strings.  Note that
 *      this is not the default hash function.  To make it the default hash
 *      function, call HashTableSetHashFunction(hashTable,
 *      HashTableStringHashFunction).
 *  ARGUMENTS:
 *      key    - the key to be hashed
 *  RETURNS:
 *      unsigned long - the unmodulated hash value of the key
\*--------------------------------------------------------------------------*/

srUInt_64 LRMHashTableStringHashFunction(const void *key) {
    const unsigned char *str = (const unsigned char *) key;
    srUInt_64 hashValue = 0;
    int i;

    for (i=0; str[i] != '\0'; i++)
        hashValue = hashValue * 37 + str[i];

    return hashValue;
}

static int LRMpointercmp(const void *pointer1, const void *pointer2) {
    return (pointer1 != pointer2);
}

static srUInt_64 LRMpointerHashFunction(const void *pointer) {
    return ((srUInt_64)( pointer - NULL )) ;
}

static int LRMisProbablePrime(srInt_64 oddNumber) {
    srInt_64 i;

    for (i=3; i<51; i+=2)
        if (oddNumber == i)
            return 1;
        else if (oddNumber%i == 0)
            return 0;

    return 1; /* maybe */
}

static srInt_64 LRMcalculateIdealNumOfBuckets(HashTable *hashTable) {
    srInt_64 idealNumOfBuckets = hashTable->numOfElements / hashTable->idealRatio;
    if (idealNumOfBuckets < 5)
        idealNumOfBuckets = 5;
    else
        idealNumOfBuckets |= 0x01; /* make it an odd number */
    while (!LRMisProbablePrime(idealNumOfBuckets))
        idealNumOfBuckets += 2;

    return idealNumOfBuckets;
}


void LRMfree_values_destroy(HashTable * tab)
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

	LRMHashTableDestroy(tab);
}

void LRMHashTableIteration(HashTable * tab, void process_item(void * key, void * hashed_obj, HashTable * tab) )
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

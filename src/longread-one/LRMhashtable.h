/*--------------------------------------------------------------------------*\
 *                   -----===== HashTable =====-----
 *
 * Author: Keith Pomakis (pomakis@pobox.com)
 * Date:   August, 1998
 * Released to the public domain.
 *
 *--------------------------------------------------------------------------
 * $Id: LRMhashtable.h,v 1.4 2019/08/28 22:31:15 cvs Exp $
\*--------------------------------------------------------------------------*/

#ifndef _HASHTABLE_H_LRM
#define _HASHTABLE_H_LRM

/* These structs should not be accessed directly from user code.
 * All access should be via the public functions declared below. */

typedef struct KeyValuePair_struct {
    const void *key;
    void *value;
    struct KeyValuePair_struct *next;
} KeyValuePair;

typedef struct {
    srInt_64 numOfBuckets;
    srInt_64 numOfElements;
    KeyValuePair **bucketArray;
    float idealRatio, lowerRehashThreshold, upperRehashThreshold;
    int (*keycmp)(const void *key1, const void *key2);
    int (*valuecmp)(const void *value1, const void *value2);
    srUInt_64 (*hashFunction)(const void *key);
    void (*keyDeallocator)(void *key);
    void (*valueDeallocator)(void *value);

    void * appendix1;
    void * appendix2;
    void * appendix3;
    srInt_64 counter1;
    srInt_64 counter2;
    srInt_64 counter3;
} HashTable;


typedef struct {
	void ** elementList;
	srInt_64 numOfElements;
	srInt_64 capacityOfElements;
	void (*elemDeallocator)(void *elem);
} ArrayList;

ArrayList * LRMArrayListCreate(int init_capacity);
void LRMArrayListDestroy(ArrayList * list);
void * LRMArrayListGet(ArrayList * list, srInt_64 n);
int LRMArrayListPush(ArrayList * list, void * new_elem);
void LRMArrayListSetDeallocationFunction(ArrayList * list,  void (*elem_deallocator)(void *elem));


void LRMHashTableIteration(HashTable * tab, void process_item(void * key, void * hashed_obj, HashTable * tab) );

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableCreate() - creates a new LRMHashTable
 *  DESCRIPTION:
 *      Creates a new LRMHashTable.  When finished with this LRMHashTable, it
 *      should be explicitly destroyed by calling the LRMHashTableDestroy()
 *      function.
 *  EFFICIENCY:
 *      O(1)
 *  ARGUMENTS:
 *      numOfBuckets - the number of buckets to start the LRMHashTable out with.
 *                     Must be greater than zero, and should be prime.
 *                     Ideally, the number of buckets should between 1/5
 *                     and 1 times the expected number of elements in the
 *                     LRMHashTable.  Values much more or less than this will
 *                     result in wasted memory or decreased performance
 *                     respectively.  The number of buckets in a LRMHashTable
 *                     can be re-calculated to an appropriate number by
 *                     calling the LRMHashTableRehash() function once the
 *                     LRMHashTable has been populated.  The number of buckets
 *                     in a LRMHashTable may also be re-calculated
 *                     automatically if the ratio of elements to buckets
 *                     passes the thresholds set by LRMHashTableSetIdealRatio().
 *  RETURNS:
 *      LRMHashTable    - a new LRMHashtable, or NULL on error
\*--------------------------------------------------------------------------*/

HashTable *LRMHashTableCreate(srInt_64 numOfBuckets);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableDestroy() - destroys an existing LRMHashTable
 *  DESCRIPTION:
 *      Destroys an existing LRMHashTable.
 *  EFFICIENCY:
 *      O(n)
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable to destroy
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableDestroy(HashTable *hashTable);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableContainsKey() - checks the existence of a key in a LRMHashTable
 *  DESCRIPTION:
 *      Determines whether or not the specified LRMHashTable contains the
 *      specified key.  Uses the comparison function specified by
 *      LRMHashTableSetKeyComparisonFunction().
 *  EFFICIENCY:
 *      O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable to search
 *      key          - the key to search for
 *  RETURNS:
 *      bool         - whether or not the specified LRMHashTable contains the
 *                     specified key.
\*--------------------------------------------------------------------------*/

int LRMHashTableContainsKey(const HashTable *hashTable, const void *key);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableContainsValue()
 *                         - checks the existence of a value in a LRMHashTable
 *  DESCRIPTION:
 *      Determines whether or not the specified LRMHashTable contains the
 *      specified value.  Unlike LRMHashTableContainsKey(), this function is
 *      not very efficient, since it has to scan linearly looking for a
 *      match.  Uses the comparison function specified by
 *      LRMHashTableSetValueComparisonFunction().
 *  EFFICIENCY:
 *      O(n)
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable to search
 *      value        - the value to search for
 *  RETURNS:
 *      bool         - whether or not the specified LRMHashTable contains the
 *                     specified value.
\*--------------------------------------------------------------------------*/

int LRMHashTableContainsValue(const HashTable *hashTable, const void *value);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTablePut() - adds a key/value pair to a LRMHashTable
 *  DESCRIPTION:
 *      Adds the specified key/value pair to the specified LRMHashTable.  If
 *      the key already exists in the LRMHashTable (determined by the comparison
 *      function specified by LRMHashTableSetKeyComparisonFunction()), its value
 *      is replaced by the new value.  May trigger an auto-rehash (see
 *      LRMHashTableSetIdealRatio()).  It is illegal to specify NULL as the
 *      key or value.
 *  EFFICIENCY:
 *      O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable to add to
 *      key          - the key to add or whose value to replace
 *      value        - the value associated with the key
 *  RETURNS:
 *      err          - 0 if successful, -1 if an error was encountered
\*--------------------------------------------------------------------------*/

int LRMHashTablePut(HashTable *hashTable, const void *key, void *value);
int LRMHashTablePutReplace(HashTable *hashTable, const void *key, void *value, int replace_key);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableGet() - retrieves the value of a key in a LRMHashTable
 *  DESCRIPTION:
 *      Retrieves the value of the specified key in the specified LRMHashTable.
 *      Uses the comparison function specified by
 *      LRMHashTableSetKeyComparisonFunction().
 *  EFFICIENCY:
 *      O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable to search
 *      key          - the key whose value is desired
 *  RETURNS:
 *      void *       - the value of the specified key, or NULL if the key
 *                     doesn't exist in the LRMHashTable
\*--------------------------------------------------------------------------*/

void * LRMHashTableGet(const HashTable *hashTable, const void *key);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableRemove() - removes a key/value pair from a LRMHashTable
 *  DESCRIPTION:
 *      Removes the key/value pair identified by the specified key from the
 *      specified LRMHashTable if the key exists in the LRMHashTable.  May trigger
 *      an auto-rehash (see LRMHashTableSetIdealRatio()).
 *  EFFICIENCY:
 *      O(1), assuming a good hash function and element-to-bucket ratio
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable to remove the key/value pair from
 *      key          - the key specifying the key/value pair to be removed
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableRemove(HashTable *hashTable, const void *key);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableRemoveAll() - removes all key/value pairs from a LRMHashTable
 *  DESCRIPTION:
 *      Removes all key/value pairs from the specified LRMHashTable.  May trigger
 *      an auto-rehash (see LRMHashTableSetIdealRatio()).
 *  EFFICIENCY:
 *      O(n)
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable to remove all key/value pairs from
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableRemoveAll(HashTable *hashTable);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableIsEmpty() - determines if a LRMHashTable is empty
 *  DESCRIPTION:
 *      Determines whether or not the specified LRMHashTable contains any
 *      key/value pairs.
 *  EFFICIENCY:
 *      O(1)
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable to check
 *  RETURNS:
 *      bool         - whether or not the specified LRMHashTable contains any
 *                     key/value pairs
\*--------------------------------------------------------------------------*/

int LRMHashTableIsEmpty(const HashTable *hashTable);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableSize() - returns the number of elements in a LRMHashTable
 *  DESCRIPTION:
 *      Returns the number of key/value pairs that are present in the
 *      specified LRMHashTable.
 *  EFFICIENCY:
 *      O(1)
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable whose size is requested
 *  RETURNS:
 *      long         - the number of key/value pairs that are present in
 *                     the specified LRMHashTable
\*--------------------------------------------------------------------------*/

srInt_64 LRMHashTableSize(const HashTable *hashTable);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableGetNumBuckets() - returns the number of buckets in a LRMHashTable
 *  DESCRIPTION:
 *      Returns the number of buckets that are in the specified LRMHashTable.
 *      This may change dynamically throughout the life of a LRMHashTable if
 *      automatic or manual rehashing is performed.
 *  EFFICIENCY:
 *      O(1)
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable whose number of buckets is requested
 *  RETURNS:
 *      long         - the number of buckets that are in the specified
 *                     LRMHashTable
\*--------------------------------------------------------------------------*/

srInt_64 LRMHashTableGetNumBuckets(const HashTable *hashTable);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableSetKeyComparisonFunction()
 *              - specifies the function used to compare keys in a LRMHashTable
 *  DESCRIPTION:
 *      Specifies the function used to compare keys in the specified
 *      LRMHashTable.  The specified function should return zero if the two
 *      keys are considered equal, and non-zero otherwise.  The default
 *      function is one that simply compares pointers.
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable whose key comparison function is being
 *                     specified
 *      keycmp       - a function which returns zero if the two arguments
 *                     passed to it are considered "equal" keys and non-zero
 *                     otherwise
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableSetKeyComparisonFunction(HashTable *hashTable,
                             int (*keycmp)(const void *key1, const void *key2));

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableSetValueComparisonFunction()
 *              - specifies the function used to compare values in a LRMHashTable
 *  DESCRIPTION:
 *      Specifies the function used to compare values in the specified
 *      LRMHashTable.  The specified function should return zero if the two
 *      values are considered equal, and non-zero otherwise.  The default
 *      function is one that simply compares pointers.
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable whose value comparison function is being
 *                     specified
 *      valuecmp     - a function which returns zero if the two arguments
 *                     passed to it are considered "equal" values and non-zero
 *                     otherwise
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableSetValueComparisonFunction(HashTable *hashTable,
                       int (*valuecmp)(const void *value1, const void *value2));

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableSetHashFunction()
 *              - specifies the hash function used by a LRMHashTable
 *  DESCRIPTION:
 *      Specifies the function used to determine the hash value for a key
 *      in the specified LRMHashTable (before modulation).  An ideal hash
 *      function is one which is easy to compute and approximates a
 *      "random" function.  The default function is one that works
 *      relatively well for pointers.  If the LRMHashTable keys are to be
 *      strings (which is probably the case), then this default function
 *      will not suffice, in which case consider using the provided
 *      LRMHashTableStringHashFunction() function.
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable whose hash function is being specified
 *      hashFunction - a function which returns an appropriate hash code
 *                     for a given key
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableSetHashFunction(HashTable *hashTable,
                              srUInt_64 (*hashFunction)(const void *key));

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableRehash() - reorganizes a LRMHashTable to be more efficient
 *  DESCRIPTION:
 *      Reorganizes a LRMHashTable to be more efficient.  If a number of
 *      buckets is specified, the LRMHashTable is rehashed to that number of
 *      buckets.  If 0 is specified, the LRMHashTable is rehashed to a number
 *      of buckets which is automatically calculated to be a prime number
 *      that achieves (as closely as possible) the ideal element-to-bucket 
 *      ratio specified by the LRMHashTableSetIdealRatio() function.
 *  EFFICIENCY:
 *      O(n)
 *  ARGUMENTS:
 *      hashTable    - the LRMHashTable to be reorganized
 *      numOfBuckets - the number of buckets to rehash the LRMHashTable to.
 *                     Should be prime.  Ideally, the number of buckets
 *                     should be between 1/5 and 1 times the expected
 *                     number of elements in the LRMHashTable.  Values much
 *                     more or less than this will result in wasted memory
 *                     or decreased performance respectively.  If 0 is
 *                     specified, an appropriate number of buckets is
 *                     automatically calculated.
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableRehash(HashTable *hashTable, srInt_64 numOfBuckets);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableSetIdealRatio()
 *              - sets the ideal element-to-bucket ratio of a LRMHashTable
 *  DESCRIPTION:
 *      Sets the ideal element-to-bucket ratio, as well as the lower and
 *      upper auto-rehash thresholds, of the specified LRMHashTable.  Note
 *      that this function doesn't actually perform a rehash.
 *
 *      The default values for these properties are 3.0, 0.0 and 15.0
 *      respectively.  This is likely fine for most situations, so there
 *      is probably no need to call this function.
 *  ARGUMENTS:
 *      hashTable    - a LRMHashTable
 *      idealRatio   - the ideal element-to-bucket ratio.  When a rehash
 *                     occurs (either manually via a call to the
 *                     LRMHashTableRehash() function or automatically due the
 *                     the triggering of one of the thresholds below), the
 *                     number of buckets in the LRMHashTable will be
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
                            float lowerRehashThreshold,
                            float upperRehashThreshold);

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableSetDeallocationFunctions()
 *              - sets the key and value deallocation functions of a LRMHashTable
 *  DESCRIPTION:
 *      Sets the key and value deallocation functions of the specified
 *      LRMHashTable.  This determines what happens to a key or a value when it
 *      is removed from the LRMHashTable.  If the deallocation function is NULL
 *      (the default if this function is never called), its reference is
 *      simply dropped and it is up to the calling program to perform the
 *      proper memory management.  If the deallocation function is non-NULL,
 *      it is called to free the memory used by the object.  E.g., for simple
 *      objects, an appropriate deallocation function may be free().
 *
 *      This affects the behaviour of the LRMHashTableDestroy(), LRMHashTablePut(),
 *      LRMHashTableRemove() and LRMHashTableRemoveAll() functions.
 *  ARGUMENTS:
 *      hashTable    - a LRMHashTable
 *      keyDeallocator
 *                   - if non-NULL, the function to be called when a key is
 *                     removed from the LRMHashTable.
 *      valueDeallocator
 *                   - if non-NULL, the function to be called when a value is
 *                     removed from the LRMHashTable.
 *  RETURNS:
 *      <nothing>
\*--------------------------------------------------------------------------*/

void LRMHashTableSetDeallocationFunctions(HashTable *hashTable,
                                       void (*keyDeallocator)(void *key),
                                       void (*valueDeallocator)(void *value));

/*--------------------------------------------------------------------------*\
 *  NAME:
 *      LRMHashTableStringHashFunction() - a good hash function for strings
 *  DESCRIPTION:
 *      A hash function that is appropriate for hashing strings.  Note that
 *      this is not the default hash function.  To make it the default hash
 *      function, call LRMHashTableSetHashFunction(HashTableStringHashFunction).
 *  ARGUMENTS:
 *      key    - the key to be hashed
 *  RETURNS:
 *      long   - the unmodulated hash value of the key
\*--------------------------------------------------------------------------*/

srUInt_64 LRMHashTableStringHashFunction(const void *key);

void LRMfree_values_destroy(HashTable * tab);
#endif /* _HASHTABLE_H */



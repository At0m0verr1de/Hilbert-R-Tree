#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
// EVERY NODE HAS BETWEEN m and M entries and children unless it is root

//--------------------------- MACRO DEFINITIONS ----------------------------//

#define M 4
#define m 2
#define MAX_CHILDREN M // M = 4; MAXIMUM NUMBER OF CHILDREN
#define MIN_CHILDREN m
#define MAX_POINTS 100
// Assumes the maximum x and y coordinates are 2^31-1
#define MAX_COORDINATE 2147483647

//--------------------------- STRUCTURE DEFINITIONS ----------------------------//

typedef struct Point Point;
struct Point
{
    int x, y;
};
typedef struct Rectangle Rectangle;
struct Rectangle
{
    Point bottom_left;
    Point top_right; // coordinates of the rectangle
    int h;           // Hilbert value of the rectangle's center
};
typedef struct LeafEntry LeafEntry;
struct LeafEntry
{
    Rectangle mbr; // Minimum bounding rectangle
    int obj_id;    // Pointer to the object description record
};
typedef struct LeafNode LeafNode;
struct LeafNode
{
    int num_entries;
    struct LeafEntry entries[MAX_CHILDREN];
};

typedef struct NonLeafEntry NonLeafEntry;
struct NonLeafEntry
{
    Rectangle mbr;             // Minimum bounding rectangle
    struct Node *child_ptr;    // Pointer to the child node
    int largest_hilbert_value; // Largest hilbert value among the data records enclosed by the MBR
};

typedef struct NonLeafNode NonLeafNode;
struct NonLeafNode
{
    int num_entries;
    struct NonLeafEntry entries[MAX_CHILDREN];
};
typedef struct Node *NODE;
struct Node
{
    int is_leaf;
    NODE parent_ptr;
    struct NonLeafNode non_leaf_node; // Non-leaf node
    struct LeafNode leaf_node;        // Leaf node
};

typedef struct HilbertRTree HilbertRTree;
struct HilbertRTree
{
    struct Node *root; // Root node of the tree
};

//--------------------------- GLOBAL VARIABLE DECLARATIONS ----------------------------//

int CURRENT_ID = 0;
int num_results = 0;
bool root_split = false;
NODE root1 = NULL;
NODE root2 = NULL;

// ---------------------------- FUNCTION DECLERATIONS ------------------------------ //
HilbertRTree *new_hilbertRTree();
NODE new_node(int is_leaf);
void InsertNode(NODE parent, NODE newNode);

uint32_t interleave(uint32_t x);
uint32_t hilbertXYToIndex(uint32_t n, uint32_t x, uint32_t y);
uint32_t hilbertValue(Rectangle rect);

NODE Insert(NODE root, Rectangle rectangle);
NODE ChooseLeaf(NODE n, Rectangle r, int h);
void AdjustTree(NODE N, NODE newNode, NODE *S, int s_size);
NODE HandleOverFlow(NODE n, Rectangle rectangle);
void preOrderTraverse(NODE n);

void printMBR(Rectangle r);
void copy_rectangle(Rectangle r1, Rectangle r2);
void adjustLHV(NODE parentNode);
void adjustMBR(NODE parentNode);
void searchGetResults(NODE root, Rectangle rectangle, LeafEntry *results);
bool isInArray(NODE *arr, int size, NODE node);
bool rectangles_equal(Rectangle *rect1, Rectangle *rect2);
bool nodes_equal(NODE node1, NODE node2);
bool intersects(Rectangle r1, Rectangle r2);
int area(Rectangle r);
int calculateLHV(NonLeafEntry entry);
int compare(const void *a, const void *b);
int find_entry_index(NODE n, Rectangle rectangle);
int numberOfSiblings(NODE *S);
Rectangle calculateMBR(Rectangle r1, Rectangle r2);
Rectangle calculateEntryMBR(NonLeafEntry entry);
NODE findLeaf(NODE root, Rectangle rectangle);
NODE *cooperatingSiblings(NODE n);
LeafEntry *search(NODE root, Rectangle rectangle);

//--------------------------- FUNCTIONS TO CREATE NEW TREE AND NODE ----------------------------//
// CREATE A NEW HILBERT R TREE STRUCTURE
HilbertRTree *new_hilbertRTree()
{
    HilbertRTree *tree = (HilbertRTree *)calloc(1,sizeof(HilbertRTree));
    tree->root = NULL;
    return tree;
}
// CREATE A NEW NODE
NODE new_node(int is_leaf)
{
    NODE node = (NODE)calloc(1,sizeof(struct Node));
    node->is_leaf = is_leaf;
    node->parent_ptr = NULL;
    node->leaf_node.num_entries = 0;
    node->non_leaf_node.num_entries = 0;
    return node;
}

// IS_LEAF = 1 IF LEAF NODE. M = MAXIMUM NUMBER OF CHILDREN THAT NODE CAN HAVE

Rectangle new_rectangle(int bottomLeft_X, int bottomLeft_y, int topRight_x, int topRight_y){
    Rectangle mbr;
    mbr.bottom_left.x = bottomLeft_X;
    mbr.bottom_left.y = bottomLeft_y;
    mbr.top_right.x = topRight_x;
    mbr.top_right.y = topRight_y;
    mbr.h = hilbertXYToIndex(5, (mbr.bottom_left.x + mbr.top_right.x) / 2, (mbr.bottom_left.y + mbr.top_right.y) / 2);
    return mbr;
}



// Function to calculate the minimum bounding rectangle that contains two given rectangles

uint32_t interleave(uint32_t x)
{
    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;
    return x;
}



// COMPUTE THE LHV OF A LEAF NODE'S PARENT NON LEAF ENTRY
int computeLeafLHV(NODE a)
{
    int LHV = 0;
    if (a == NULL)
    {
        return LHV;
    }
    for (int i = 0; i < a->leaf_node.num_entries; i++)
    {
        if (a->leaf_node.entries[i].mbr.h > LHV)
        {
            LHV = a->leaf_node.entries[i].mbr.h;
        }
    }
    return LHV;
}




// COMPARE TWO NON LEAF ENTRIES BASED ON LHV
int compareNonLeafEntry(const void *a, const void *b)
{
    const struct NonLeafEntry *s1 = a;
    const struct NonLeafEntry *s2 = b;
    return s1->largest_hilbert_value - s2->largest_hilbert_value;
}

// compare: Rectangle
int compare(const void *a, const void *b)
{
    const struct LeafEntry *s1 = a;
    const struct LeafEntry *s2 = b;
    return s1->mbr.h - s2->mbr.h;
}


// PRE-ORDER TRAVERSAL OF A NODE; PRINTING A INTERNAL NODE AND THEN
void preOrderTraverse(NODE n)
{
    // IF ALL NODES TRAVERSED
    if (!n)
    {
        return;
    }
    // IF NODE IS A LEAF: PRINT ALL ENTRIES
    if (n->is_leaf == 1)
    {

        for (int i = 0; i < n->leaf_node.num_entries; i++)
        {
            printf("Leaf Node Entry: %d\n", i); //
            printf("Object_ID = %d: ", n->leaf_node.entries[i].obj_id);
            printMBR(n->leaf_node.entries[i].mbr);
        }
    }
    // IF NODE IS NON-LEAF
    else
    {
        for (int i = 0; i < n->non_leaf_node.num_entries; i++)
        {
            printf("Internal node Entry %d\n", i);
            printMBR(n->non_leaf_node.entries[i].mbr);
            preOrderTraverse(n->non_leaf_node.entries[i].child_ptr);
        }
    }
}


uint32_t hilbertXYToIndex(uint32_t n, uint32_t x, uint32_t y)
{
    x = x << (16 - n);
    y = y << (16 - n);

    uint32_t A, B, C, D;

    // Initial prefix scan round, prime with x and y
    {
        uint32_t a = x ^ y;
        uint32_t b = 0xFFFF ^ a;
        uint32_t c = 0xFFFF ^ (x | y);
        uint32_t d = x & (y ^ 0xFFFF);

        A = a | (b >> 1);
        B = (a >> 1) ^ a;

        C = ((c >> 1) ^ (b & (d >> 1))) ^ c;
        D = ((a & (c >> 1)) ^ (d >> 1)) ^ d;
    }

    {
        uint32_t a = A;
        uint32_t b = B;
        uint32_t c = C;
        uint32_t d = D;

        A = ((a & (a >> 2)) ^ (b & (b >> 2)));
        B = ((a & (b >> 2)) ^ (b & ((a ^ b) >> 2)));

        C ^= ((a & (c >> 2)) ^ (b & (d >> 2)));
        D ^= ((b & (c >> 2)) ^ ((a ^ b) & (d >> 2)));
    }

    {
        uint32_t a = A;
        uint32_t b = B;
        uint32_t c = C;
        uint32_t d = D;

        A = ((a & (a >> 4)) ^ (b & (b >> 4)));
        B = ((a & (b >> 4)) ^ (b & ((a ^ b) >> 4)));

        C ^= ((a & (c >> 4)) ^ (b & (d >> 4)));
        D ^= ((b & (c >> 4)) ^ ((a ^ b) & (d >> 4)));
    }

    // Final round and projection
    {
        uint32_t a = A;
        uint32_t b = B;
        uint32_t c = C;
        uint32_t d = D;

        C ^= ((a & (c >> 8)) ^ (b & (d >> 8)));
        D ^= ((b & (c >> 8)) ^ ((a ^ b) & (d >> 8)));
    }

    // Undo transformation prefix scan
    uint32_t a = C ^ (C >> 1);
    uint32_t b = D ^ (D >> 1);

    // Recover index bits
    uint32_t i0 = x ^ y;
    uint32_t i1 = b | (0xFFFF ^ (i0 | a));

    return ((interleave(i1) << 1) | interleave(i0)) >> (32 - 2 * n);
}

// RETURNS TRUE IF ALL NODES IN S HAVE M NUMBER OF ENTRIES; FALSE IF ATLEAST ONE OF THE NDOES IN S HAS LESS THAN M ENTRIES
bool allNodesFull(NODE *S, int numSiblings)
{
    bool allFull = true;
    for (int i = 0; i < numSiblings; i++)
    {
        if (S[i]->is_leaf == 1)
        {
            if (S[i]->leaf_node.num_entries < M)
            {
                allFull = false;
                break;
            }
        }
        else if (S[i]->is_leaf == 0)
        {
            if (S[i]->non_leaf_node.num_entries < M)
            {
                allFull = false;
                break;
            }
        }
    }
    return allFull;
}

 
// CALCULATE THE MBR RECTANGLES FROM TWO RECTANGLES: POTENTIAL HELPER FUNCTION
Rectangle calculateMBR(Rectangle r1, Rectangle r2)
{
    Point bottom_leftNew;
    Point top_rightNew;
    if (r1.bottom_left.x <= r2.bottom_left.x)
    {
        bottom_leftNew.x = r1.bottom_left.x;
    }
    else
    {
        bottom_leftNew.x = r2.bottom_left.x;
    }
    if (r1.bottom_left.y <= r2.bottom_left.y)
    {
        bottom_leftNew.y = r1.bottom_left.y;
    }
    else
    {
        bottom_leftNew.y = r2.bottom_left.y;
    }

    if (r1.top_right.x <= r2.top_right.x)
    {
        top_rightNew.x = r2.top_right.x;
    }
    else
    {
        top_rightNew.x = r1.top_right.x;
    }
    if (r1.top_right.y <= r2.top_right.y)
    {
        top_rightNew.y = r2.top_right.y;
    }
    else
    {
        top_rightNew.y = r1.top_right.y;
    }

    Rectangle new_rect;
    new_rect.bottom_left = bottom_leftNew;
    new_rect.top_right = top_rightNew;
    return new_rect;
}

// CALCULATE THE MBR FOR A NON LEAF ENTRY BASED ON ENTRIES IN ITS CHILD NODE
Rectangle calculateEntryMBR(NonLeafEntry entry)
{
    // FOR EACH NON LEAF ENTRY; CALCULATE MBR FROM CHILD  NODES
    // FIND -> LOWEST X, LOWEST Y AND HIGHEST X, HIGHEST Y
    Rectangle mbr;
    NODE next_node = entry.child_ptr;
    int low_x = 0;
    int low_y = 0;
    int high_x = 0;
    int high_y = 0;

    // IF CHILD NODE IS A LEAF; GET MBR FOR ALL ITS ENTRIES
    if (next_node != NULL && next_node->is_leaf == 1)
    {
        high_x = next_node->leaf_node.entries[0].mbr.top_right.x;
        high_y = next_node->leaf_node.entries[0].mbr.top_right.y;
        low_x = next_node->leaf_node.entries[0].mbr.top_right.x;
        low_y = next_node->leaf_node.entries[0].mbr.top_right.y;

        for (int i = 0; i < next_node->leaf_node.num_entries; i++)
        {
            Rectangle obj_mbr = next_node->leaf_node.entries[i].mbr;
            low_x = (obj_mbr.bottom_left.x <= low_x) ? obj_mbr.bottom_left.x : low_x;
            low_y = (obj_mbr.bottom_left.y <= low_y) ? obj_mbr.bottom_left.y : low_y;
            high_x = (obj_mbr.top_right.x >= high_x) ? obj_mbr.top_right.x : high_x;
            high_y = (obj_mbr.top_right.y >= high_y) ? obj_mbr.top_right.y : high_y;
        }
    }
    // IF CHILD NODE IS NON-LEAF;
    else if (next_node != NULL && next_node->is_leaf == 0)
    {
        // NON LEAF NODE: ASSUMING IT HAS BEEN RECURSED UPON IN ADJUSTMBR; CALCULATE MBR FROM ITS ENTRIES
        low_x = next_node->non_leaf_node.entries[0].mbr.top_right.x;
        low_y = next_node->non_leaf_node.entries[0].mbr.top_right.y;
        high_x = next_node->non_leaf_node.entries[0].mbr.top_right.x;
        high_y = next_node->non_leaf_node.entries[0].mbr.top_right.y;

        for (int i = 0; i < next_node->non_leaf_node.num_entries; i++)
        {
            Rectangle child_mbr = next_node->non_leaf_node.entries[i].mbr;
            low_x = (child_mbr.bottom_left.x <= low_x) ? child_mbr.bottom_left.x : low_x;
            low_y = (child_mbr.bottom_left.y <= low_y) ? child_mbr.bottom_left.y : low_y;
            high_x = (child_mbr.top_right.x >= high_x) ? child_mbr.top_right.x : high_x;
            high_y = (child_mbr.top_right.y >= high_y) ? child_mbr.top_right.y : high_y;
        }
    }
    // SET THE COORDINATES FOR THE MBR
    mbr = new_rectangle(low_x, low_y, high_x, high_y);
    return mbr;
}

// PRINT THE MBR - TOP RIGHT POINT AND BOTTOM LEFT POINT: POTENTIAL HELPER FUNCTION
void printMBR(Rectangle rect)
{

    printf("MBR = (%d, %d) to  (%d, %d)\n", rect.bottom_left.x, rect.bottom_left.y, rect.top_right.x, rect.top_right.y);
    return;
}

// STORE ALL ELAF ENTRIES OF LEAF NODES IN S INTO E
void store_all_leaf_entries(LeafEntry *E, int *num_entries, NODE *S, int numSiblings)
{
    //FOR EVERY SIBLING NODE
    for (int i = 0; i < numSiblings; i++)
    {
        // SIBLING NODE IS A LEAF
        if (S[i]->is_leaf == 1)
        {
            // STORES THE MBR OF ENTRIES INTO THE E ARRAY
            for (int j = 0; j < S[i]->leaf_node.num_entries; j++)
            {
                E[(*num_entries)++] = S[i]->leaf_node.entries[j];
            }
        }
    }
}


// FUNCTION TO CREATE A NEW LEAF ENTRY WITH THE RECTANGLE PASSED AND A NEW OBJECT ID
LeafEntry new_leafentry(Rectangle rectangle)
{
    LeafEntry le;
    le.mbr = rectangle;
    le.obj_id = ++CURRENT_ID;
    return le;
}

// CREATE A NEW NON LEAF-ENTRY GIVEN THE CHILD NODE
NonLeafEntry new_nonleafentry(NODE newNode)
{

    NonLeafEntry nle;
    nle.child_ptr = newNode;
    // IMP: SET THE MBR AND LHV OF THE NEW ENTRY AFTER DEFINING CHILD POINTER
    nle.mbr = calculateEntryMBR(nle);
    nle.largest_hilbert_value = calculateLHV(nle);
    return nle;
}

// STORE ALL NON-LEAF ENTRIES IN THE NODES IN S INTO SET E (FOR HANDLEOVERFLOWNODE)
void store_all_nonleaf_entries(NonLeafEntry *E, NODE *S, int *num_entries, int numSiblings)
{
    // FOR EACH SIBLING NODE
    for (int i = 0; i < numSiblings; i++)
    {
        // SIBLING NODE IS NON-LEAF
        if (S[i]->is_leaf == 0)
        {
            // STORE ALL THE NON LEAF ENTRIES OF THE NODE
            for (int j = 0; j < S[i]->non_leaf_node.num_entries; j++)
            {
                E[(*num_entries)++] = S[i]->non_leaf_node.entries[j];
            }
        }
    }
}

void distribute_nonleaf_entries_evenly(NODE *S, int numSiblings, NonLeafEntry *E, int *num_entries)
{
    int num_entries_per_node = (*num_entries) / numSiblings;
    int remainder_entries = (*num_entries) % numSiblings;

    // FOR EACH SIBLING NODE

    // DISTRIBUTION LIST[I] = NUMBER OF ENTRIES FOR THE ITH SIBLINGS
    int *distributionList = (int *)calloc(numSiblings, sizeof(int));

    //EACH SIBLING HAS ATLEAST NUM_ENTRIES_PER_NODE ENTRIES
    for (int i = 0; i < numSiblings; i++)
    {
        distributionList[i] = num_entries_per_node;
    }
    //ADD THE REMAINDER ENTRIES
    for (int j = 0; j < remainder_entries; j++)
    {
        distributionList[j]++;
    }


    int done = 0;
    //ADD THE ENTRIES TO THE SIBLINGS
    for (int j = 0; j < numSiblings; j++)
    {
        //SET THE NON LEAF ENTRIES TO 0
        S[j]->non_leaf_node.num_entries = 0;

        for (int l = 0; l < distributionList[j]; l++)
        {

            //ADD THE NON LEAF ENTRY
            S[j]->non_leaf_node.entries[l] = E[done];
            E[done].child_ptr->parent_ptr = S[j];
            S[j]->non_leaf_node.num_entries++;
            done++;
        }
        //SET REMAINING ENTRIES TO 0
        for (int l = distributionList[j]; l < M; l++)
        {
            S[j]->non_leaf_node.entries[l].mbr.bottom_left.x = 0;
            S[j]->non_leaf_node.entries[l].mbr.bottom_left.y = 0;
            S[j]->non_leaf_node.entries[l].mbr.top_right.x = 0;
            S[j]->non_leaf_node.entries[l].mbr.top_right.y = 0;
            S[j]->non_leaf_node.entries[l].mbr.h = 0;
            S[j]->non_leaf_node.entries[l].child_ptr = NULL;
            S[j]->non_leaf_node.entries[l].largest_hilbert_value = 0;
        }
    }

    free(distributionList);
    return;
}

NODE HandleOverFlowNode(NODE parentNode, NODE new_node1)
{

    printf("HANDLE OVERFLOW NODE CALLED\n");

    // TO INSERT NEW_NODE AS A CHILD POINTER IN A NON LEAF ENTRY
    NonLeafEntry entry = new_nonleafentry(new_node1);

    // SET OF COOPERATING SIBLINGS FOR THE PARENTNODE
    NODE *S = cooperatingSiblings(parentNode);
 
    // NUMBER OF SIBLINGS
    int numSiblings = numberOfSiblings(S);

    //SET THE POSSIBLE NODE UPON SLITTING OF PARENTNODE TO NULL
    NODE NN = NULL;

    // SET OF ALL NON LEAF ENTRIES IN PARENTNODE AND SIBLINGS
    int *num_entries = (int *)malloc(sizeof(int));

    *num_entries = 0;
    NonLeafEntry *E = (NonLeafEntry *)calloc(MAX_POINTS, sizeof(NonLeafEntry));

    // ADD ALL ENTRIES TO E FROM THE SIBLINGS
    store_all_nonleaf_entries(E, S, num_entries, numSiblings);

    // ADD NEW NON LEAF ENTRY TO E
    E[(*num_entries)++] = entry;

    // SORT THE SET OF NON LEAF ENTRIES BASED ON LHV OF NON LEAF ENTRIES
    qsort(E, *num_entries, sizeof(NonLeafEntry), compareNonLeafEntry);

    // IF ALL SIBLINGS ARE FULL OR NOT
    bool allFull = true;
    allFull = allNodesFull(S, numSiblings); // True if all nodes in S are full

    // IF ALL SIBLINGS ARE FULL
    if (allFull)
    {
        printf("ALL PARENT NODE'S %d SIBLINGS ARE FULL\n", numSiblings);

        // CREATE A NEW NODE
        NN = new_node(0);

        if (parentNode->parent_ptr == NULL)
        {
            // PARENT NODE IS ROOT: ROOT WAS SPLIT
            root_split = true;
        }

        // ADD NN TO SIBLINGS
        S[numSiblings++] = NN; // ADD THE NEW NODE TO THE SET OF SIBLINGS

        // printf("PARENT NODE: NUMBER OF NON LEAF ENTRIES = %d NUMBER OF SIBLINGS = %d\n", *num_entries, numSiblings);

    }

    distribute_nonleaf_entries_evenly(S, numSiblings, E, num_entries);
    free(E); free(num_entries);
    return NN;
    
}

// HELPER FUNCTION TO RETURN TRUE IF TWO RECTANGLES HAVE THE SAME BOTTOM LEFT AND TOP RIGHT POINTS ARE SAME
bool rectangles_equal(Rectangle *rect1, Rectangle *rect2)
{
    if (rect1->bottom_left.x != rect2->bottom_left.x ||
        rect1->bottom_left.y != rect2->bottom_left.y ||
        rect1->top_right.x != rect2->top_right.x ||
        rect1->top_right.y != rect2->top_right.y)
    {
        return false;
    }
    return true;
}

// HELPER FUNCTION TO CHECK IF TWO NODES ARE THE SAME; LEAF/NONLEAF AND IDENTICAL ENTRIS
bool nodes_equal(NODE node1, NODE node2)
{
    if (node1->is_leaf != node2->is_leaf)
    {
        return false;
    }

    if (node1->is_leaf)
    {
        // Compare LeafNode fields
        struct LeafNode *leaf1 = &node1->leaf_node;
        struct LeafNode *leaf2 = &node2->leaf_node;
        if (leaf1->num_entries != leaf2->num_entries)
        {
            return false;
        }
        for (int i = 0; i < leaf1->num_entries; i++)
        {
            if (rectangles_equal(&leaf1->entries[i].mbr, &leaf2->entries[i].mbr) &&
                leaf1->entries[i].obj_id == leaf2->entries[i].obj_id)
            {
                continue;
            }
            else
            {
                return false;
            }
        }
    }
    else
    {
        // Compare NonLeafNode fields
        struct NonLeafNode *non_leaf1 = &node1->non_leaf_node;
        struct NonLeafNode *non_leaf2 = &node2->non_leaf_node;
        if (non_leaf1->num_entries != non_leaf2->num_entries)
        {
            return false;
        }
        for (int i = 0; i < non_leaf1->num_entries; i++)
        {
            if (rectangles_equal(&non_leaf1->entries[i].mbr, &non_leaf2->entries[i].mbr) &&
                non_leaf1->entries[i].child_ptr == non_leaf2->entries[i].child_ptr &&
                non_leaf1->entries[i].largest_hilbert_value == non_leaf2->entries[i].largest_hilbert_value)
            {
                continue;
            }
            else
            {
                return false;
            }
        }
    }
    return true;
}

// CALCULATE LHV FOR A NON LEAF ENTRY: BY EVAUATING ALL CHILD PTR ETRIES
int calculateLHV(NonLeafEntry entry)
{
    int max_h = 0;
    NODE node = entry.child_ptr;
    // CALCULATE MAXIMUM H OF NODE
    if (node->is_leaf == 1)
    {
        max_h = node->leaf_node.entries[0].mbr.h;
        for (int i = 0; i < node->leaf_node.num_entries; i++)
        {
            if (node->leaf_node.entries[i].mbr.h > max_h)
            {
                max_h = node->leaf_node.entries[i].mbr.h;
            }
        }
    }
    else if (node->is_leaf == 0)
    {
        max_h = node->non_leaf_node.entries[0].largest_hilbert_value;
        // NON LEAF CHILD NODE
        for (int i = 0; i < node->non_leaf_node.num_entries; i++)
        {
            if (node->non_leaf_node.entries[i].largest_hilbert_value > max_h)
            {
                max_h = node->non_leaf_node.entries[i].largest_hilbert_value;
            }
        }
    }
    return max_h;
}


// HELPER FUNCTION TO CHECK IN ARRAY
bool isInArray(NODE *arr, int size, NODE node)
{
    for (int i = 0; i < size; i++)
    {
        if (nodes_equal(arr[i], node))
        {
            return true;
        }
    }
    return false;
}



/*ADJUST TREE ASCEND FROM LEAF TOWARDS ROOT AND ADJUST MBR AND LHV VALUES*/
void adjustLHV(NODE parentNode)
{
    if (parentNode == NULL)
    {
        return;
    }

    for (int i = 0; i < parentNode->non_leaf_node.num_entries; i++)
    {
        parentNode->non_leaf_node.entries[i].largest_hilbert_value = calculateLHV(parentNode->non_leaf_node.entries[i]);
    }
    adjustLHV(parentNode->parent_ptr);
}



void adjustMBR(NODE parentNode)
{
    if (parentNode == NULL)
    {
        return;
    }

    for (int i = 0; i < parentNode->non_leaf_node.num_entries; i++)
    {
        parentNode->non_leaf_node.entries[i].mbr = calculateEntryMBR(parentNode->non_leaf_node.entries[i]);
    }
    adjustMBR(parentNode->parent_ptr);
}


// HELPER FUNCTION TO INSERT THE NEWNODE AS THE CHILD NODE OF A NONLEAF ENTRY INTO NODE PARENT
void InsertNode(NODE parent, NODE newNode)
{
    // SET THE PARENT POINTER OF ADDED NODE
    newNode->parent_ptr = parent;

    // NON LEAF ENTRY CREATED WITH CHILD NODE AS NEWNODE
    NonLeafEntry entry = new_nonleafentry(newNode);

    // INSERT THE NEW NON LEAF ENTRY ACCORDING TO LHV
    int i = 0;
    while (i < parent->non_leaf_node.num_entries && parent->non_leaf_node.entries[i].largest_hilbert_value <= entry.largest_hilbert_value)
    {
        i++;
    }

    // SHIFT VALUES IN PARENT NODE
    for (int j = parent->non_leaf_node.num_entries; j > i; j--)
    {
        parent->non_leaf_node.entries[j] = parent->non_leaf_node.entries[j - 1];
    }

    // INSERT ACCORDING TO HILBERT VALUE
    parent->non_leaf_node.entries[i] = entry;
    parent->non_leaf_node.num_entries++;
}

// FUNCTION TO FIND A LEAF WHERE RECTANGLE R SHOULD BE INSERTED
NODE ChooseLeaf(NODE n, Rectangle r, int h) // PARAMETERS: NODE N, RECTANGLE R, HILBERT VALUE FO CENTRE OF RECTANGLE: H
{
    /* RETURNS THE LEAF NODE IN WHICH TO PLACE A NEW RECTANGE*/
    /* IF N IS A LEAF, RETURN N*/
    // printf("CHOOSE LEAF CALLED\n");

    //C2. LEAF CHECK: IF N IS A LEAF NODE; RETURN N
    if (n->is_leaf == 1)
    {
        return n;
    }

    //INITIALISE VARIABLES TO STORE CURRENT NODE SELECTED AND MINIMUM LHV (GREATER THAN H OF RECTANGLE)
    float min_LHV = n->non_leaf_node.entries[0].largest_hilbert_value;
    NODE next_node = NULL;

    // C3. CHOOSE SUBTREE: IF N IS A NON-LEAF NODE, CHOOSE ENTRY (R, PTR, LHV) WITH MINIMUM LHV GREATER THAN H
    for (int i = 0; i < n->non_leaf_node.num_entries; i++)
    {
        if (n->non_leaf_node.entries[i].largest_hilbert_value >= h && n->non_leaf_node.entries[i].largest_hilbert_value <= min_LHV)
        {
            min_LHV = n->non_leaf_node.entries[i].largest_hilbert_value;
            next_node = n->non_leaf_node.entries[i].child_ptr;
        }
    }


    /* IF ALL CHILDREN HAVE LHV LESS THAN H */
    //NOTE: NOW MIN_LHV STORES THE LARGEST LHV ENCOUNTERED YET AS ALL ENTRIES HAVE LHV LESS THAN H
    if (next_node == NULL)
    {
        min_LHV = n->non_leaf_node.entries[0].largest_hilbert_value;
        // CHOOSE THE CHILD NODE WITH LARGEST LHV (AS ALL HAVE LHV LESS THAN HILBERT VALUE OF RECTANGLE)
        for (int i = 0; i < n->non_leaf_node.num_entries; i++)
        {
            if (n->non_leaf_node.entries[i].largest_hilbert_value >= min_LHV)
            {
                min_LHV = n->non_leaf_node.entries[i].largest_hilbert_value;
                next_node = n->non_leaf_node.entries[i].child_ptr;
            }
        }
    }


    // C4. DESCEND UNTIL A LEAF NODE IS REACHED
    return ChooseLeaf(next_node, r, h);

}


// SORT THE SIBLING NODES IN S USING BUBBLE SORT
void sortSiblings(NODE *S, int numSiblings)
{
    for (int i = 0; i < numSiblings - 1; i++)
    {
        for (int j = 0; j < numSiblings - i - 1; j++)
        {
            // SORTING HAPPENS BASED ON THE LARGEST HILBERT VALUES OF THEIR ENTRIES
                if (computeLeafLHV(S[j]) > computeLeafLHV(S[j + 1]))
                {
                    NODE temp = S[j];
                    S[j] = S[j + 1];
                    S[j + 1] = temp;
                }
            

        }
    }
}


// HELPER FUNCTION TO RETURN AN ARRAY OF COOPERATING SIBLING NODES AND THE NODE ITSELF (AS REQD IN MULTIPLE ALGORITHMS): S=2; 1 COOPERATING SIBLINGS
NODE *cooperatingSiblings(NODE n)
{
    NODE *S = (NODE *)calloc(MAX_POINTS, sizeof(NODE));

    // INITIALISE ALL ENTIRES TO NULL
    for (int i = 0; i < MAX_POINTS; i++)
    {
        S[i] = NULL;
    }

    // ADD THE NODE ITSELF TO THE SET
    S[0] = n;

    // SET INDEX TO 0
    int numSiblingsCP = 0;

    // PARENTNODE IS THE PARENT NODE OF N
    NODE parentNode = n->parent_ptr;

    // IF PARENTNODE IS NULL;  N IS ROOT: NO 1 SIBLINGS
    if (parentNode == NULL)
    {
        return S;
    }

    // GO FROM CURRENT NODE TO ROOT FINDING SIBLINGS

    // 1. FIND INDEX OF THE NODE IN THE PARENT NODE
    int index = -1;
    for (int i = 0; i < parentNode->non_leaf_node.num_entries; i++)
    {
        if (parentNode->non_leaf_node.entries[i].child_ptr == n)
        {
            index = i;
            break;
        }
    }

    // IF NODE ON LEFT IS AVAILABLE
    if (index > 0)
    {
        S[0] = parentNode->non_leaf_node.entries[index - 1].child_ptr;
        S[1] = n;
        numSiblingsCP++;
    }

    // IF NODE ON RIGHT IS AVAILABLE
    if (index < parentNode->non_leaf_node.num_entries - 1)
    {
        S[++numSiblingsCP] = parentNode->non_leaf_node.entries[index + 1].child_ptr;
    }

    // SORT S ON COMPUTELHV(ELEMENT): LHV OF THE NON LEAF ENTRIES
    sortSiblings(S, numSiblingsCP + 1);

    return S;
}

// FUNCTION TO FIND THE NUMBER OF NODES IN ARRAY S; USED TO FIND THE NUMBER OF SIBLINGS
int numberOfSiblings(NODE *S)
{

    int numSiblings = 0;

    for (int i = 0; i < MAX_POINTS; i++)
    {

        if (S[i] != NULL)
        {

            numSiblings++;

        }
    }

    return numSiblings;
}

// HELPER FUNCTION TO COPY BOTTOM LEFT POINT AND TOP RIGHT POINT OF RECTANGLE R2 TO RECTANGLE R1
void copy_rectangle(Rectangle r1, Rectangle r2)
{
    r1.bottom_left.x = r2.bottom_left.x;
    r1.top_right.x = r2.top_right.x;
    r1.bottom_left.y = r2.bottom_left.y;
    r1.top_right.y = r2.top_right.y;
}


// HELPER FUNCTION TO DISTRIBUTE THE LEAF ENTRIES IN ARRAY E INTO THE NODES IN ARRAY S EVENLY
void distribute_leaf_entries_evenly(NODE *S, int numSiblings, LeafEntry *E, int *num_entries)
{
    int num_entries_per_node = (*num_entries) / numSiblings;
    int remainder_entries = (*num_entries) % numSiblings;

    // DISTRIBUTION LIST[I] = NUMBER OF LEAF ENTRIES FOR THE ITH SIBLING LEAF NODE
    int *distributionList = (int *)calloc(numSiblings, sizeof(int));

    for (int i = 0; i < numSiblings; i++)
    {
        distributionList[i] = num_entries_per_node;
    }
    for (int j = 0; j < remainder_entries; j++)
    {
        distributionList[j]++;
    }

    // DISTRIBUTE THE LEAF ENTRIES AMONGST THE SIBLINGS
    int done = 0;


    // FOR EACH SIBLINGS
    for (int j = 0; j < numSiblings; j++)
    {
        // SET NUMBER OF LEAF ENTRIES TO 0
        S[j]->leaf_node.num_entries = 0;
        for (int l = 0; l < distributionList[j]; l++)
        {

            S[j]->leaf_node.entries[l] = E[done];
            S[j]->leaf_node.num_entries++;
            done++;
        }
        for (int l = distributionList[j]; l < M; l++)
        {
            S[j]->leaf_node.entries[l].mbr.bottom_left.x = 0;
            S[j]->leaf_node.entries[l].mbr.bottom_left.y = 0;
            S[j]->leaf_node.entries[l].mbr.top_right.x = 0;
            S[j]->leaf_node.entries[l].mbr.top_right.y = 0;
            S[j]->leaf_node.entries[l].mbr.h = 0;
            S[j]->leaf_node.entries[l].obj_id = 0;
        }
    }
    free(distributionList);
    return;
}

// TRYING TO RETAIN OBJECT_IDS BY STORING LEAF ENTRIES INSTEAD OF RECTANGLES
NODE HandleOverFlow(NODE n, Rectangle rectangle)
{

    NODE *S = cooperatingSiblings(n);      // CONTAINS COOPERATING SIBLINGS AND NODE

    int numSiblings = numberOfSiblings(S); // SIZE OF SET S

    int *num_entries = (int *)malloc(sizeof(int));
    *num_entries = 0; // SIZE OF SET E

    NODE NN = NULL;

    // H1: SET OF ALL ENTRIES FROM SIBLINGS
    // E = SET OF ALL ENTRIES FROM N AND S-1 COOPERATING SIBLINGS [ALL STORED IN SET S]
    LeafEntry *E = (LeafEntry *)calloc(MAX_POINTS, sizeof(LeafEntry));
    store_all_leaf_entries(E, num_entries, S, numSiblings);

    // CREATE NEW LEAF ENTRY WITH RECTANGLE
    LeafEntry entry = new_leafentry(rectangle);

    // H2. ADD R TO E; ADD NEW LEAF ENTRY TO SET E
    E[(*num_entries)++] = entry;

    // SORT THE ENTRIES IN E BASED ON THEIR HILBERT VALUE
    qsort(E, *(num_entries), sizeof(LeafEntry), compare);

    // CHECK IF ANY NODE IS NOT FULL
    bool allFull = true;
    allFull = allNodesFull(S, numSiblings);

    // H4. IF ATLEAST ONE OF SIBLINGS IS NOT FULL
    if (allFull)
    {
        // DISTRIBUTE E EVENLY AMONG THE S+1 NODES ACCORDING TO THE HILBERT VALUE
        // CREATE A NEW NODE NN

        // printf("ALL SIBLINGS OF LEAF NODE ARE FULL$d\n");

        NN = new_node(1); // CREATE A NEW LEAF NODE;

        // IF ROOT WAS SPLIT TO CREATE NEW NODE
        if (n->parent_ptr == NULL)
        {
            // ROOT NODE IS SPLIT
            root_split = true;
        }

        // ADD NN TO SIBLINGS
        S[numSiblings++] = NN; // ADD THE NEW NODE TO THE SET OF SIBLINGS

    }
    
    
    distribute_leaf_entries_evenly(S, numSiblings, E, num_entries);
    free(E);
    free(num_entries);
    return NN;

}

// RETURNS TRUE IF ROOT NODE WAS SPLIT
void AdjustTree(NODE N, NODE NN, NODE *S, int s_size)
{
    // STOP IF ROOT LEVEL REACHED
    NODE Np = N->parent_ptr;
    NODE new_node = NULL;


    // printf("ADJUST TREE CALLED\n");

    // PARENT = NULL; ROOT LEVEL
    //A1. IF ROOT LEVEL IS REACHED; STOP
    if (Np == NULL)
    {
        // printf("ROOT LEVEL REACHED IN ADJUST TREE\n");
        return;
    }

    // NN == NULL; SPLITTING DID NOT HAPPEN
    if (NN == NULL)
    {
        // printf("NO NEW NODE PASSED TO ADJUST TREE\n");
    }

    // INSERT SPLIT NODE INTO PARENT
    //A2. PROPOGATE NODE SPLIT UPWARD
    if (NN != NULL)
    {

        // printf("HANDLING NEW NODE IN ADJUST TREE\n");

        // INSERT IN CORRECT ORDER IF ROOM IN PARENT NODE
        if (Np->non_leaf_node.num_entries < MAX_CHILDREN)
        {

            // printf("PARENT NODE HAS SPACE\n");
            InsertNode(Np, NN);

        }
        //OTHERWISE INNVOKE HANDLEOVERFLOWNODE
        else
        {
            // PARENT NODE MUST BE SPLIT
            // printf("PARENT NODE IS FULL\n");

            // HANDLEOVERFLOWNODE: WHEN PARENT NODE MUST BE SPLIT
            //IF NP IS SPLIT; NEW_NODE IS THE NEW NODE
            new_node = HandleOverFlowNode(Np, NN);

            if (new_node != NULL)
            {

                // printf("PARENT NODE WAS SPLIT\n");

            }

            // IF ROOT NODE WAS SPLIT BY HANDLEOVERFLOW
            if (Np->parent_ptr == NULL && new_node != NULL)
            {

                // printf("ROOT SHOULD HAVE SPLIT!\n");
                root_split = true;
                root1 = Np;
                root2 = new_node;

            }
        }
    }

    // A3. ADJUST MBR AND LHV IN PARENT LEVEL
    // P = SET OF PARENT NODES FOR NODES IN S

    NODE *P = (NODE *)calloc(MAX_POINTS, sizeof(NODE));

    int numParents = 0;

    for (int i = 0; i < s_size; i++)
    {

        if (S[i]->parent_ptr == NULL)
        {

            continue;

        }

        NODE parent = S[i]->parent_ptr;

        if (parent)
        {

            P[numParents++] = parent;

        }

    }

    // ADJUST MBR AND LHV VALUES
    for (int i = 0; i < numParents; i++)
    {

        //A3. ADJUST CORRESPONDINNG MBRs AND LHVs OF NODES IN P
        NODE parent = P[i];
        adjustMBR(parent);
        adjustLHV(parent);

    }

    // A4. MOVE UP TO NEXT LEVEL: NN = NEW_NODE, S = P, NUMSIBLINGS = NUMPARENTS
    AdjustTree(Np, new_node, P, numParents);

}

// INSERTION: INSERT A RECTANGLE INTO THE RTREE BY PASSING ITS ROOT NODE AND THE RECTANGLE
NODE Insert(NODE root, Rectangle rectangle)
{
    //I1. GET THE SUITABLE LEAF NODE WHERE RECTANGLE SHOULD BE INSERTED
    NODE leafNode = ChooseLeaf(root, rectangle, rectangle.h);
    //SET THE POSSIBLE NODE UPON SPLITTING TO NULL
    NODE newLeafNode = NULL;
    //SET OF COOPERATING SIBLINGS AND NUMBER OF SIBLING NDOES
    NODE *S = cooperatingSiblings(leafNode);
    int numSiblings = numberOfSiblings(S);


    //I2. INSERT R IN A LEAF NODE L: IF L HAS AN EMPTY SLOT
    if (leafNode->leaf_node.num_entries < M)
    {
        // LEAF NODE HAS SPACE
        // printf("LEAF NODE HAS SPACE: %d\n", rectangle.h);


        //FIND INDEX WHERE RECTAGNGLE SHOULD BE INSERTED AS PER HILBERT VALUE
        int i = 0;
        while (i < leafNode->leaf_node.num_entries && leafNode->leaf_node.entries[i].mbr.h < rectangle.h)
        {
            i++;
        }

        //SHIFT ENTRIES TO ADD THE NEW LEAF ENTRY
        for (int j = leafNode->leaf_node.num_entries; j > i; j--)
        {
            leafNode->leaf_node.entries[j] = leafNode->leaf_node.entries[j - 1];
        }


        // CREATE THE LEAF ENTRY WITH NEW RECTANGLE
        LeafEntry entry = new_leafentry(rectangle);

        // INSERT THE LEAF ENTRY INTO NODE
        leafNode->leaf_node.entries[i] = entry;
        leafNode->leaf_node.num_entries++;

    }

    //I2. INSERT R IN A LEAF NODE L; IF L IS FULL
    else if (leafNode->leaf_node.num_entries >= M)
    {
        // LEAF NODE IS FULL
        // printf("LEAF NODE IS FULL: %d\n", rectangle.h);

        //INVOKE HANDLEOVERFLOW; RETURNING A NEW NODE IS SPLIT HAPPENED
        newLeafNode = HandleOverFlow(leafNode, rectangle);

        //IF HANDLEOVERFLOW RETURNED A NODE
        if (newLeafNode)
        {

            // printf("LEAF NODE SPLIT: %d\n", rectangle.h);
            //NEW NODE IS A LEAF NODE
            newLeafNode->is_leaf = 1;
            //I3. SET S SHOULD CONTAIN L, COOPERATING SIBLINGS AND THE NEW NODE 
            S[numSiblings++] = newLeafNode;

        }
    }

    //I3. PROPOGATE CHANGES UPWARD: INVOKE ADJUSTTREE
    root_split = false;
    AdjustTree(leafNode, newLeafNode, S, numSiblings);
 
    //I4. GROW TREE TALLER: IF ROOT WAS SPLIT (AND ROOT WAS THE ONLY NODE -> LEAF NODE)
    if (leafNode->parent_ptr == NULL && newLeafNode != NULL)
    {

            // printf("NEW ROOT FORMED");
            root_split = true;
            root1 = leafNode;
            root2 = newLeafNode;

    }
    // RETURNS THE NEW LEAF IF SPLIT WAS INEVITABLE
    
    //I4. GROW TREE TALLER; IF ADJUSTTREE AND NODE SPLIT PROPOGATION CAUSED THE ROOT TO SPLIT
    if (root_split)
    {

        // printf("ROOT IS SPLIT: %d\n", rectangle.h);
        //FORM NEW ROOT WHICH IS NOT A LEAF NODE
        NODE newRoot = new_node(0);
        newRoot->non_leaf_node.num_entries = 2;

        //CREATE THE NEW NON LEAF ENTRIES; WITH ROOT1 AND ROOT2 AS CHILD NODES
        NonLeafEntry entry1 = new_nonleafentry(root1);
        NonLeafEntry entry2 = new_nonleafentry(root2);

        //ADD THE NON LEAF ENTRIES TO THE NEWLY CREATED ROOT
        newRoot->non_leaf_node.entries[0] = entry1;
        newRoot->non_leaf_node.entries[1] = entry2;
        
        //SET THE PARENT POINTERS OF ROOT1 AND ROOT2
        root1->parent_ptr = newRoot;
        root2->parent_ptr = newRoot;

        free(S);
        return newRoot;
    }
    free(S);
    return root;
}

// SEARCH ALGORITHM
//-> NONLEAF - THOSE WITH MBR INTERSECTING THE QUERY WINDOW W
//-> LEAF - THOSE WITH MBR INTERSECTING THE QUERY WINDOW W
/*ALL RECTANGLES THAT OVERLAP A SEARCH RECTANGLE*/
bool intersects(Rectangle r1, Rectangle r2)
{
    return !(
        r1.top_right.x < r2.bottom_left.x || r2.top_right.x < r1.bottom_left.x ||
        r1.top_right.y < r2.bottom_left.y || r2.top_right.y < r1.bottom_left.y);
}



// SEARCH SHOULD RETURN AN ARRAY OF RESULTS
void searchGetResults(NODE root, Rectangle rectangle, LeafEntry *results)
{
    if (root->is_leaf == 1)
    {
        for (int i = 0; i < root->leaf_node.num_entries; i++)
        {
            if (intersects(root->leaf_node.entries[i].mbr, rectangle))
            {
                results[num_results++] = root->leaf_node.entries[i];
            }
        }
    }
    else
    {
        for (int i = 0; i < root->non_leaf_node.num_entries; i++)
        {
            if (intersects(root->non_leaf_node.entries[i].mbr, rectangle))
            {
                searchGetResults(root->non_leaf_node.entries[i].child_ptr, rectangle, results);
            }
        }
    }
}



LeafEntry *search(NODE root, Rectangle rectangle)
{
    // NUMBER OF RESULTS: STORED IN GLOBAL VARIABLE
    num_results = 0;
    LeafEntry *results = (LeafEntry *)calloc(MAX_POINTS, sizeof(LeafEntry));
    searchGetResults(root, rectangle, results);
    return results;
}
// handle overflow



/*FIND THE LEAF NODE CONTAINING A RECTANGLE R*/
NODE findLeaf(NODE root, Rectangle rectangle)
{
    if (root->is_leaf == 1)
    {
        return root;
    }

    /* IF NON LEAF; FIND OVERLAPPING ENTRY*/
    int i;
    for (i = 0; i < root->non_leaf_node.num_entries; i++)
    {
        if (intersects(root->non_leaf_node.entries[i].mbr, rectangle))
        {
            break;
        }
    }
    return findLeaf(root->non_leaf_node.entries[i].child_ptr, rectangle);
}



int find_entry_index(NODE n, Rectangle rectangle)
{
    int index = -1;
    for (int i = 0; i < n->leaf_node.num_entries; i++)
    {
        if (n->leaf_node.entries[i].mbr.bottom_left.x == rectangle.bottom_left.x &&
            n->leaf_node.entries[i].mbr.bottom_left.y == rectangle.bottom_left.y &&
            n->leaf_node.entries[i].mbr.top_right.x == rectangle.top_right.x &&
            n->leaf_node.entries[i].mbr.top_right.y == rectangle.top_right.y)
        {
            index = i;
            break;
        }
    }
    return index;
}



void print_mbr(Rectangle r)
{
    printf("Top right point: (%d, %d)\n", r.top_right.x, r.top_right.y);
    printf("Bottom left point: (%d, %d)\n", r.bottom_left.x, r.bottom_left.y);
}



int main()
{
    FILE *fp;
    Point points[MAX_POINTS];
    int num_points = 0;

    //OPENING FILE CONTAINING THE INPUT DATA POINTS
    fp = fopen("data1.txt", "r");
    

    //IF FILE COULD NOT BE OPENED
    if (fp == NULL)
    {
        printf("Error opening file.\n");
        return 1;
    }


    // INPUT DATA POINTS FROM THE FILE: NUMBER OF POINTS SPECIFIED BY MACRO MAX_POINTS
    for (int i = 0; i < MAX_POINTS; i++)
    {
        fscanf(fp, "%d %d\n", &points[num_points].x, &points[num_points].y);
        // printf("POINT  = %d %d\n", points[num_points].x, points[num_points].y);
        ++num_points;
    }

    //CLOSING THE FILE AFTER READING INPUT
    fclose(fp);

    // LETS ASSUME THE RECTANGLES INITIALLY ARE THE DATA POINTS WITH THE BOTTOM LEFT CORNER = TOP RIGHT CORNER = DATA-POINT
    Rectangle rectangles[MAX_POINTS];


    // CREATE RECTANGLES FROM THE DATA POINTS
    for (int i = 0; i < MAX_POINTS; i++)
    {
        rectangles[i] = new_rectangle(points[i].x, points[i].y, points[i].x, points[i].y);
    }

    // QUICK SORT RECTANGLES ARRAY
    qsort(rectangles, MAX_POINTS, sizeof(struct Rectangle), compare); // ARGUMENTS = ARRAY, NUMBER OF ELEMENTS, SIZE OF EACH ELEMENT, COMPARISON FUNCTION

    // PRINT THE H-VALUES OF RECTANGLES
    // for (int i = 0; i < MAX_POINTS; i++)
    // {
    //     printf("RECTANGLE WITH H-VALUE = %d\n", rectangles[i].h);
    // }

    // CREATE A HILBERT R TREE
    HilbertRTree *Rtree = new_hilbertRTree();

    // INITIALLY THE ROOT OF TREE IS A LEAF;
    Rtree->root = new_node(1);
 
    // THE PARENT POINTER OF ROOT IS NULL
    Rtree->root->parent_ptr = NULL;

    // INSERT THE RECTANGLES INTO THE HILBERT RTREE
    for (int i = 0; i < MAX_POINTS; i++)
    {

        Rtree->root = Insert(Rtree->root, rectangles[i]);
    }

    // PERFORM A PRE ORDER TRAVERSAL OF THE TREE
    preOrderTraverse(Rtree->root);
    return 0;
}

#include <stdint.h>
#include <stdio.h>

// EVERY NODE HAS BETWEEN m and M entries and children unless it is root
#define M 4
#define m 2
#define MAX_CHILDREN M // M = 4; MAXIMUM NUMBER OF CHILDREN
#define MIN_CHILDREN m
#define MAX_POINTS 5
#define MAX_COORDINATE 2147483647  // Assumes the maximum x and y coordinates are 2^31-1

typedef struct Point Point;
typedef struct Rectangle Rectangle;
typedef struct LeafEntry LeafEntry;
typedef struct NonLeafEntry NonLeafEntry;
typedef struct Node* NODE;
typedef struct HilbertRTree HilbertRTree;

int CURRENT_ID = 0;
int num_results = 0;
bool root_split = false;
int numberOfElements = 0;

NODE root1 = NULL;
NODE root2 = NULL;




struct Point{
    int x,y;
};

struct Rectangle {
    Point bottom_left;
    Point top_right;    // Coordinates of the rectangle
    float h;    // Hilbert value of the rectangle center
};

struct LeafEntry {
    Rectangle mbr;  // Minimum bounding rectangle
    int obj_id;  // Pointer to the object description record
};

struct LeafNode {
    int num_entries;
    struct LeafEntry entries[MAX_CHILDREN];
};

struct NonLeafEntry {
    Rectangle mbr;  // Minimum bounding rectangle
    NODE child_ptr;  // Pointer to the child node
    int largest_hilbert_value;  // Largest hilbert value among the data records enclosed by the MBR
};

struct NonLeafNode {
    int num_entries;
    struct NonLeafEntry entries[MAX_CHILDREN];    
};

struct Node {
    int is_leaf;
    NODE parent_ptr;
    union {
        struct NonLeafNode non_leaf_node;  // Non-leaf node
        struct LeafNode leaf_node;  // Leaf node
    } u;
};

struct HilbertRTree {
    int height;  // Height of the tree
    NODE root;  // Root node of the tree
};

// ----------------------------FUNCTION DECLERATIONS------------------------------
HilbertRTree* new_hilbertRTree();
NODE new_node(int is_leaf);
void InsertNode(NODE parent, NODE newNode);

uint32_t interleave(uint32_t x);
uint32_t hilbertXYToIndex(uint32_t n, uint32_t x, uint32_t y);
uint32_t hilbertValue(Rectangle rect);

NODE Insert(NODE root, Rectangle rectangle);
void delete(NODE root, Rectangle rectangle);
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
Rectangle calculateMBR(Rectangle r1, Rectangle r2);
Rectangle calculateEntryMBR(NonLeafEntry entry);
NODE findLeaf(NODE root, Rectangle rectangle);
NODE *cooperatingSiblings(NODE n, int *numSiblings);
LeafEntry *search(NODE root, Rectangle rectangle);

// -----------------------FUNCTIONS TO CREATE NODE AND TREE-----------------------

HilbertRTree* new_hilbertRTree(){
    HilbertRTree* tree = (HilbertRTree*)malloc(sizeof(HilbertRTree));
    tree->root = NULL;
    tree->height = 0;
    return tree;
}

NODE new_node(int is_leaf){
    NODE node = (NODE)malloc(sizeof(struct Node));
    node->is_leaf = is_leaf;
    node->parent_ptr = NULL;
    if(is_leaf){
        node->u.leaf_node.num_entries = 0;
    }
    else{
        node->u.non_leaf_node.num_entries = 0;
    }
    return node;
}

void InsertNode(NODE parent, NODE newNode){
    // newNode->is_leaf=0;
    // parent->is_leaf=0;
    newNode->parent_ptr = parent;
    parent->non_leaf_node.num_entries++;
    parent->non_leaf_node.entries[(parent->non_leaf_node.num_entries) - 1].child_ptr = newNode;
    parent->non_leaf_node.entries[(parent->non_leaf_node.num_entries) - 1].mbr = calculateEntryMBR(parent->non_leaf_node.entries[(parent->non_leaf_node.num_entries) - 1]);
    parent->non_leaf_node.entries[(parent->non_leaf_node.num_entries) - 1].largest_hilbert_value = calculateLHV(parent->non_leaf_node.entries[(parent->non_leaf_node.num_entries) - 1]);
}


// -----------------------Functions to calculate Hilbert value-----------------------

// Helper function to calculate the Hilbert value of a point
uint32_t interleave(uint32_t x){
    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;
    return x;
}

// Calculate the Hilbert value of a point
uint32_t hilbertXYToIndex(uint32_t n, uint32_t x, uint32_t y){
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

// Calculate the Hilbert value of a rectangle
uint32_t hilbertValue(Rectangle rect){
    uint32_t x = (rect.bottom_left.x + rect.top_right.x)/2;
    uint32_t y = (rect.bottom_left.y + rect.top_right.y)/2;
    uint32_t hilbert_value = hilbertXYToIndex(16, x, y);
    return hilbert_value;
}

// -----------------------Functions to operate on the Hilbert Tree-------------------

// Inserting a rectangle into the tree
NODE Insert(NODE root, Rectangle rectangle){
    NODE leafNode = ChooseLeaf(root, rectangle, rectangle.h);
    NODE newLeafNode = NULL;

    if (leafNode->leaf_node.num_entries < M){   // LEAF NODE HAS SPACE
        int i = 0;

        // While the hilbert value of the current entry is less than the hilbert value of the rectangle
        // And the current entry is not the last entry
        while ((i < leafNode->leaf_node.num_entries) && (leafNode->leaf_node.entries[i].mbr.h < rectangle.h)){
            i++;
        }

        // Shift all the entries to the right by 1
        for (int j = leafNode->leaf_node.num_entries; j > i; j--){ 
            leafNode->leaf_node.entries[j] = leafNode->leaf_node.entries[j - 1];
        }

        // INSERT ACCORDING TO HILBERT ORDER AND RETURN
        LeafEntry entry;
        entry.mbr = rectangle;
        entry.obj_id = ++CURRENT_ID;
        leafNode->leaf_node.entries[i] = entry;
        leafNode->leaf_node.num_entries++;

        root_split = false;
        int* numSiblings = malloc(sizeof(int));
        NODE *S = cooperatingSiblings(leafNode, numSiblings);
        AdjustTree(leafNode, newLeafNode, S, ++numSiblings);

        free(numSiblings);

        return root;
    }
    else
    {
        // LEAF NODE IS FULL
        newLeafNode = HandleOverFlow(leafNode, rectangle);
        newLeafNode->is_leaf = 1;
        if (leafNode->parent_ptr == NULL)
        {
            int* numSiblings = malloc(sizeof(int));
            NODE *S = cooperatingSiblings(leafNode, numSiblings);
            AdjustTree(leafNode, newLeafNode, S, numSiblings + 1);

            NODE newRoot = (NODE)malloc(sizeof(struct Node));
            newRoot->is_leaf = 0;
            newRoot->non_leaf_node.num_entries = 2;
            newRoot->parent_ptr = NULL;

            newRoot->non_leaf_node.entries[0].child_ptr = leafNode;
            newRoot->non_leaf_node.entries[0].largest_hilbert_value = calculateLHV(newRoot->non_leaf_node.entries[0]);
            newRoot->non_leaf_node.entries[0].mbr = calculateEntryMBR(newRoot->non_leaf_node.entries[0]);

            newRoot->non_leaf_node.entries[1].child_ptr = newLeafNode;
            newRoot->non_leaf_node.entries[1].largest_hilbert_value = calculateLHV(newRoot->non_leaf_node.entries[1]);
            newRoot->non_leaf_node.entries[1].mbr = calculateEntryMBR(newRoot->non_leaf_node.entries[1]);

            leafNode->parent_ptr = newRoot;
            newLeafNode->parent_ptr = newRoot;
            return newRoot;
        }
        // RETURNS THE NEW LEAF IF SPLIT WAS INEVITABLE
    }

    // Propogate changes upward
    // FORM A SET S CONTAINING L: COOPERATING SIBLINGS AND NEW LEAF (IF ANY)

    root_split = false;
    int* numSiblings = malloc(sizeof(int));
    NODE *S = cooperatingSiblings(leafNode, numSiblings);
    AdjustTree(leafNode, newLeafNode, S, numSiblings + 1);

    // IF NODE SPLIT CAUSED ROOT TO SPLIT, CREATEA NEWROOT WITH CHILDREN
    // AS RESULTING NODES
    //  Check if the root split
    if (root_split)
    {
        NODE newRoot = (NODE)malloc(sizeof(struct Node));
        newRoot->is_leaf = 0;
        newRoot->non_leaf_node.num_entries = 2;
        newRoot->parent_ptr = NULL;

        newRoot->non_leaf_node.entries[0].child_ptr = root1;
        newRoot->non_leaf_node.entries[0].largest_hilbert_value = calculateLHV(newRoot->non_leaf_node.entries[0]);
        newRoot->non_leaf_node.entries[0].mbr = calculateEntryMBR(newRoot->non_leaf_node.entries[0]);

        newRoot->non_leaf_node.entries[1].child_ptr = root2;
        newRoot->non_leaf_node.entries[1].largest_hilbert_value = calculateLHV(newRoot->non_leaf_node.entries[1]);
        newRoot->non_leaf_node.entries[1].mbr = calculateEntryMBR(newRoot->non_leaf_node.entries[1]);

        root1->parent_ptr = newRoot;
        root2->parent_ptr = newRoot;
        return newRoot;
    }
    return root;
}

/*DELETE(RECTANGLE R)*/
void delete(NODE root, Rectangle rectangle){
    // D1. FIND THE HOST LEAF
    // FIND THE HOST LEAF: N = LEAF CONTAINING ENTRY WITH RECTANGLE
    NODE n = findLeaf(root, rectangle);
    NODE parentNode = n->parent_ptr;
    NODE *S = (NODE *)malloc(sizeof(NODE) * MAX_POINTS);
    int* numSiblings = malloc(sizeof(int));

    // D2. DELETE R: REMOVE R FROM NODE N
    // FIND INDEX OF ENTRY
    int index = find_entry_index(n, rectangle);

    if (index == -1)
    {
        printf("ENTRY NOT FOUND TO DELETE");
        return;
    }

    // ENTRY FOUND; INDEX>=0 DELETE ENTRY AT INDEX I
    for (int i = index; i < n->leaf_node.num_entries - 1; i++)
    {
        n->leaf_node.entries[i] = n->leaf_node.entries[i + 1];
    }
    n->leaf_node.num_entries--;

    // D3. IF NODE UNDERFLOWS: LESS THAN m ENTRIES
    if (n->leaf_node.num_entries < MIN_CHILDREN)
    {
        // BORROW ENTRIES FROM COOPERATING SIBLINGS
        S = cooperatingSiblings(n, numSiblings);
        numSiblings++;

        int num_borrowed = 0;
        bool allUnderflow = false;
        int i = 1;

        // ARE ALL COOPERATING SIBLINGS READY TO UNDERFLOW?
        for (int i = 1; i < numSiblings; i++)
        {
            if (S[i]->leaf_node.num_entries > MIN_CHILDREN)
            {
                allUnderflow = true;
            }
        }
        // IF NOT ALL READY TO UNDERFLOW: BORROW FROM S COOPERATING SIBLINGS
        if (!allUnderflow)
        {
            // WHILE NODE IS IN UNDERFLOW AND ALL SIBLINGS NOT EXPLORED
            while (n->leaf_node.num_entries < MIN_CHILDREN && i < numSiblings)
            {
                // IF SIBLING IS NOT NULL AND HAS MORE THAN MINIMUM CHILDREN
                while (S[i] != NULL && S[i]->leaf_node.num_entries > MIN_CHILDREN && n->leaf_node.num_entries < MIN_CHILDREN)
                {
                    // LOGIC TO TRANSFER/BORROW AN ENTRY: BORROW 1ST ENTRY
                    n->leaf_node.entries[n->leaf_node.num_entries].mbr = S[i]->leaf_node.entries[0].mbr;
                    n->leaf_node.entries[n->leaf_node.num_entries].obj_id = S[i]->leaf_node.entries[0].obj_id;
                    n->leaf_node.num_entries++;
                    // SHIFT ALL ENTRIES ONE POSITION TO LEFT
                    for (int j = 0; j < S[i]->leaf_node.num_entries - 1; j++)
                    {
                        S[i]->leaf_node.entries[j] = S[i]->leaf_node.entries[j + 1];
                    }
                    S[i]->leaf_node.num_entries--;
                    num_borrowed++;
                    break;
                }
                i++;
            }
        }
        // IF ALL SIBLINGS ARE READY TO UNDERFLOW
        if (allUnderflow)
        {
            // MERGE THE CURRENT NODE WITH A SIBLING NODE
            for (int i = 0; i < n->leaf_node.num_entries; i++)
            {
                // INSERT THE ENTRIES TO A COOPERATING SIBLING
                S[1]->leaf_node.entries[S[i]->leaf_node.num_entries] = n->leaf_node.entries[i];
                S[1]->leaf_node.num_entries++;
            }
            parentNode = n->parent_ptr;
            free(n);

            // DELETE THE ENTRY CONTAINING CHILD_PTR TO DELETED NODE
            for (int i = index; i < parentNode->non_leaf_node.num_entries - 1; i++)
            {
                parentNode->non_leaf_node.entries[i] = parentNode->non_leaf_node.entries[i + 1];
            }
            parentNode->non_leaf_node.num_entries--;
            // UPDATE LHV AND MBR OF PARENT NODE AND ENTRY
            adjustLHV(parentNode);
            adjustMBR(parentNode);
        }
        // BORROW SOME ENTRIES FROM S COOPERATING SIBLINGS Having more than minimum
        // IF ALL SIBLINGS READS TO UNDERFLOW; MERGE S+1 TO S NODES
        // ADJUST THE RESULTING NODES
    }
    adjustLHV(parentNode);
    adjustMBR(parentNode);
    // ADJUST MBR AND LHV IN PARENT LEVELS
    // FORM A SET S CONTAINING L AND COOPERATING SIBLINGS [IF UNDERFLOW HAD OCCURED]
    AdjustTree(n, NULL, S, numSiblings);
    free(S);
}

// Chooses a lead for inserting the rectangle r with hilbert value h
NODE ChooseLeaf(NODE n, Rectangle r, int h){
    /* RETURNS THE LEAF NODE IN WHICH TO PLACE A NEW RECTANGE*/
    /* IF N IS A LEAF, RETURN N*/
    if (n->is_leaf == 1)
    {
        return n;
    }
    /* IF N IS A NON-LEAF NODE, CHOOSE ENTRY (R, PTR, LHV) WITH MINIMUM LHV GREATER THAN H*/
    float min_LHV = INFINITY;
    NODE next_node = NULL;
    for (int i = 0; i < n->non_leaf_node.num_entries; i++)
    {
        if (n->non_leaf_node.entries[i].largest_hilbert_value > h && n->non_leaf_node.entries[i].largest_hilbert_value < min_LHV)
        {
            min_LHV = n->non_leaf_node.entries[i].largest_hilbert_value;
            next_node = n->non_leaf_node.entries[i].child_ptr;
        }
    }
    /* IF ALL CHILDREN HAVE LHV LESS THAN H */
    if (next_node == NULL)
    {
        // CHOOSE THE CHILD NODE WITH LARGEST LHV
        for (int i = 0; i < n->non_leaf_node.num_entries; i++)
        {

            if (n->non_leaf_node.entries[i].largest_hilbert_value > min_LHV)
            {
                min_LHV = n->non_leaf_node.entries[i].largest_hilbert_value;
                next_node = n->non_leaf_node.entries[i].child_ptr;
            }
        }
    }

    /* DESCEND UNTIL A LEAF NODE IS REACHED*/
    return ChooseLeaf(next_node, r, h);
}

// Adjusts the tree after insertion of the new rectangle
// Returns true if root node was split
void AdjustTree(NODE N, NODE newNode, NODE *S, int s_size){

    NODE parentNode = N->parent_ptr;
    NODE new_node = NULL;

    // Stop if root level reached
    if (!parentNode){
        return;
    }

    // Insert split node in parent node
    if (newNode){
        if (parentNode->non_leaf_node.num_entries < MAX_CHILDREN){ // Parent Node has space
            InsertNode(parentNode, newNode);
        } else {
            // PARENT NODE MUST BE SPLIT
            // new_node = HandleOverFlowNode(Np, NN);
            // IF ROOT NODE WAS SPLIT BH HANDLEOVERFLOW
            if (parentNode->parent_ptr == NULL && new_node != NULL)
            {
                root_split = true;
                root1 = parentNode;
                root2 = new_node;
            }
        }
    }

    // ADJUST MBR AND LHV IN PARENT LEVEL
    // P = SET OF PARENT NODES FOR NODES IN S
    NODE *P = (NODE *)malloc(sizeof(struct Node) * MAX_POINTS);
    int numParents = 0;
    for (int i = 0; i < s_size; i++){
        NODE parent = S[i]->parent_ptr;
        if (parent && !isInArray(P, numParents, parent)){

            P[numParents++] = parent;
        }
    }
    // ADJUST MBR AND LHV VALUES
    for (int i = 0; i < numParents; i++){
        NODE parent = P[i];
        adjustMBR(parent);
        adjustLHV(parent);
    }

    // NEXT LEVEL
    AdjustTree(parentNode, new_node, P, numParents);
}

// Handle overflow of a node on insertion of a new rectangle in the node
NODE HandleOverFlow(NODE n, Rectangle rectangle){
    // E = SET OF ALL ENTRIES FROM N AND S-1 COOPERATING SIBLINGS

    int* numSiblings = malloc(sizeof(int));    // SIZE OF SET S
    NODE *S = (NODE *)malloc(sizeof(NODE) * MAX_POINTS);

    S = cooperatingSiblings(n, numSiblings);   // CONTAINS COOPERATING SIBLINGS AND NODE
    numSiblings++;

    int num_entries = 0;                       // SIZE OF SET E
    Rectangle E[MAX_POINTS];

    // H1: SET OF ALL ENTRIES FROM SIBLINGS
    for (int i = 0; i < numSiblings; i++){
        // SIBLING NODE IS A LEAF
        if (S[i]->is_leaf == 1){
            // STORES THE MBR OF ENTRIES INTO THE E ARRAY
            for (int j = 0; j < M; j++){
                E[num_entries++] = S[i]->leaf_node.entries[j].mbr;
            }
        }
        // SIBLING NODE IS NON-LEAF
        else{
            // STORES THE MBR OF ENTRIES INTO THE E ARRAY
            for (int j = 0; j < M; j++){
                E[num_entries++] = S[i]->non_leaf_node.entries[j].mbr;
            }
        }
    }
    // ADD R TO E
    E[num_entries++] = rectangle;

    // CHECK IF ANY NODE IS NOT FULL
    bool allFull = true;
    for (int i = 0; i < numSiblings; i++){
        if (S[i]->is_leaf == 1){
            if (S[i]->leaf_node.num_entries < M){
                allFull = false;
            }
        }
        else{
            if (S[i]->non_leaf_node.num_entries < M){
                allFull = false;
            }
        }
    }

    // IF ATLEAST ONE OF SIBLINGS IS NOT FULL
    if (!allFull){
        /*num_entries: total number of entries*/
        /*num_siblings+1: total number of siblings to distribute in */

        // SORT THE ENTRIES IN E BASED ON THEIR HILBERT VALUE
        int n = sizeof(E) / sizeof(E[0]);
        qsort(E, n, sizeof(struct Rectangle), compare);

        int num_entries_per_node = num_entries / numSiblings;
        int remainder_entries = num_entries % numSiblings;

        // FOR EACH SIBLING NODE
        int currIndex = 0;
        int k = 0;
        for (int i = 0; i < numSiblings; i++){
            int startIndex = 0;
            S[i]->leaf_node.num_entries = 0;
            for (k = currIndex; k < currIndex + num_entries_per_node; k++){
                // INSERT THE ENTRY
                S[i]->leaf_node.entries[startIndex++].mbr = E[k];
                S[i]->leaf_node.entries[startIndex].obj_id = ++CURRENT_ID;
                S[i]->leaf_node.num_entries++;
            }
        }
        
        // ADDING REMAINDER ENTRIES
        int j = 0;

        for (int i = 0; i < remainder_entries; i++){
            if (S[j]->leaf_node.num_entries < MAX_CHILDREN){
                S[j]->leaf_node.entries[S[j]->leaf_node.num_entries].mbr = E[k++];
                S[j]->leaf_node.entries[S[j]->leaf_node.num_entries].obj_id = ++CURRENT_ID;
                S[j]->leaf_node.num_entries++;
            }
            else{
                j++;
            }
        }
        // DISTRIBUTE E EVENLY AMONG THE S NODES ACCORDING TO THE HILBERT VALUE
    }
    // IF ALL COOPERATING SIBLINGS ARE FULL
    else {
        // DISTRIBUTE E EVENLY AMONG THE S+1 NODES ACCORDING TO THE HILBERT VALUE
        // CREATE A NEW NODE NN
        NODE NN = (NODE)malloc(sizeof(struct Node));
        if (n->parent_ptr == NULL){
            root_split = true;
        }
        S[numSiblings] = NN; // ADD THE NEW NODE TO THE SET OF SIBLINGS
        numSiblings++;

        int n = sizeof(E) / sizeof(E[0]);
        // qsort(E, n, sizeof(struct Rectangle), compare);
        int num_entries_per_node = num_entries / numSiblings;
        int remainder_entries = num_entries % numSiblings;
        // FOR EACH SIBLING NODE
        int *distributionList = (int *)calloc(numSiblings, sizeof(int));

        for (int i = 0; i < numSiblings; i++){
            distributionList[i] = num_entries_per_node;
        }
        for (int j = 0; j < remainder_entries; j++){
            distributionList[j]++;
        }
        int a = distributionList[0];
        int b = distributionList[1];
        int currIndex = 0;
        int k = 0;

        int done = 0;
        for (int j = 0; j < numSiblings; j++)
        {
            int startIndex = 0;
            S[j]->leaf_node.num_entries = 0;
            for (int l = 0; l < distributionList[j]; l++){
                S[j]->leaf_node.entries[l].mbr = E[done];
                // printf("Updated j='%d'and done='%d' and l='%d'\n",j,done,l);
                S[j]->leaf_node.entries[l].obj_id = ++CURRENT_ID;
                S[j]->leaf_node.num_entries++;
                done++;
            }
            for (int l = distributionList[j]; l < M; l++){
                S[j]->leaf_node.entries[l].mbr.bottom_left.x = 0;
                S[j]->leaf_node.entries[l].mbr.bottom_left.y = 0;
                S[j]->leaf_node.entries[l].mbr.top_right.x = 0;
                S[j]->leaf_node.entries[l].mbr.top_right.y = 0;
                S[j]->leaf_node.entries[l].mbr.h = 0;
                S[j]->leaf_node.entries[l].obj_id = 0;
                // printf("Went through j='%d' and l='%d'\n",j,l);
            }
        }

        // ADDING REMAINDER ENTRIES
        // int j = 0;
        // for (int i = 0; i < remainder_entries; i++)
        // {
        //     if (S[j]->leaf_node.num_entries < MAX_CHILDREN)
        //     {
        //             S[j]->leaf_node.entries[S[j]->leaf_node.num_entries].mbr = E[k++];
        //             S[j]->leaf_node.entries[S[j]->leaf_node.num_entries].obj_id = ++CURRENT_ID;
        //             S[j]->leaf_node.num_entries++;
        //     }
        //     else
        //     {
        //             j++;
        //     }
        // }

        return NN;
    }
    return NULL;
}

// Pre Order Traversal of the tree and print the nodes
void preOrderTraverse(NODE n) {
    if (n->is_leaf == 1)
    {
        printf("Leaf node\n");
        for (int i = 0; i < n->leaf_node.num_entries; i++)
        {
            printf("Object_ID = %d: ", n->leaf_node.entries[i].obj_id);
            printMBR(n->leaf_node.entries[i].mbr);
        }
    }
    else
    {
        printf("Internal node\n");
        for (int i = 0; i < n->non_leaf_node.num_entries; i++)
        {
            printMBR(n->non_leaf_node.entries[i].mbr);
            traverse(n->non_leaf_node.entries[i].child_ptr);
        }
    }
}


// ------------------------ HELPER FUNCTIONS ------------------------

// Helper Function: Print MBR
void printMBR(Rectangle r) {
    printf("Bottom Left: (%d, %d), Top Right: (%d, %d)\n", r.bottom_left.x, r.bottom_left.y, r.top_right.x, r.top_right.y);
}

// Helper Function: Copy Rectangle r2 into r1
void copy_rectangle(Rectangle r1, Rectangle r2) {
    r1.bottom_left.x = r2.bottom_left.x;
    r1.top_right.x = r2.top_right.x;
    r1.bottom_left.y = r2.bottom_left.y;
    r1.top_right.y = r2.top_right.y;
}

// Helper Function: Finds the index of the entry in the leaf node
int find_entry_index(NODE n, Rectangle rectangle) {
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

// Helper Function: Checks if node is present in the given array of nodes 
bool isInArray(NODE *arr, int size, NODE node) {
    for (int i = 0; i < size; i++)
    {
        if (nodes_equal(arr[i], node))
        {
            return true;
        }
    }
    return false;
}

// Helper Function: Checks if rectangles r1 and r2 are the same rectangles
bool rectangles_equal(Rectangle *rect1, Rectangle *rect2) {
    if (rect1->bottom_left.x != rect2->bottom_left.x ||
        rect1->bottom_left.y != rect2->bottom_left.y ||
        rect1->top_right.x != rect2->top_right.x ||
        rect1->top_right.y != rect2->top_right.y)
    {
        return false;
    }
    return true;
}

// Helper Function: Checks if two nodes are equal
bool nodes_equal(NODE node1, NODE node2) {
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

// Helper Function: Calculates the area of the rectangle
int area(Rectangle r) {
    return (r.top_right.x - r.bottom_left.x) * (r.top_right.y - r.bottom_left.y);
}

// Helper Function: Calculates the minimum bounding rectangle that contains two given rectangles
Rectangle calculateMBR(Rectangle r1, Rectangle r2) {
    Point bottom_leftNew;
    Point top_rightNew;
    if(r1.bottom_left.x <= r2.bottom_left.x){
        bottom_leftNew.x = r1.bottom_left.x;
    }
    else{
        bottom_leftNew.x = r2.bottom_left.x;
    }
    if(r1.bottom_left.y <= r2.bottom_left.y){
        bottom_leftNew.y = r1.bottom_left.y;
    }
    else{
        bottom_leftNew.y = r2.bottom_left.y;
    }

    if(r1.top_right.x <= r2.top_right.x){
        top_rightNew.x = r2.top_right.x;
    }
    else{
        top_rightNew.x = r1.top_right.x;
    }
    if(r1.top_right.y <= r2.top_right.y){
        top_rightNew.y = r2.top_right.y;
    }
    else{
        top_rightNew.y = r1.top_right.y;
    }

    Rectangle new_rect;
    new_rect.bottom_left = bottom_leftNew;
    new_rect.top_right = top_rightNew;
    return new_rect;
}

// Helper Function: Calcuates the LHV of the MBR for both leaf and non-leaf nodes
int calculateLHV(NonLeafEntry entry) {
    int max_h;
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
        return max_h;
    }
    else
    {
        max_h = node->non_leaf_node.entries[0].mbr.h;
        // NON LEAF CHILD NODE
        for (int i = 0; i < node->non_leaf_node.num_entries; i++)
        {
            if (node->non_leaf_node.entries[i].mbr.h > max_h)
            {
                max_h = node->non_leaf_node.entries[i].mbr.h;
            }
        }
        return max_h;
    }
}

// Helper Function: returns the difference between two hilbert values
int compare(const void *a, const void *b) {
    const struct Rectangle *s1 = a;
    const struct Rectangle *s2 = b;
    return s1->h - s2->h;
}

// Helper Function: Searches all rectangles in the tree that intersect with the query window
bool intersects(Rectangle r1, Rectangle r2) {
    return !(
        r1.top_right.x < r2.bottom_left.x || r2.top_right.x < r1.bottom_left.x ||
        r1.top_right.y < r2.bottom_left.y || r2.top_right.y < r1.bottom_left.y);
}

// Helper Function: Finds the node containing the rectangle
NODE findLeaf(NODE root, Rectangle rectangle) {
    if (root->is_leaf == 1){
        return root;
    }

    /* IF NON LEAF; FIND OVERLAPPING ENTRY*/
    int i;
    for (i = 0; i < root->non_leaf_node.num_entries; i++){
        if (intersects(root->non_leaf_node.entries[i].mbr, rectangle)){
            break;
        }
    }
    return findLeaf(root->non_leaf_node.entries[i].child_ptr, rectangle);
}

// Helper Function: Uses searchGetResults() to return a LeafEntry type array of 
//                      all rectangles that intersect with the query window
LeafEntry *search(NODE root, Rectangle rectangle) {
    // NUMBER OF RESULTS: STORED IN GLOBAL VARIABLE
    num_results = 0;
    LeafEntry *results = (LeafEntry *)malloc(sizeof(LeafEntry) * MAX_POINTS);
    searchGetResults(root, rectangle, results);
    return results;
}


// Helper Function:
void searchGetResults(NODE root, Rectangle rectangle, LeafEntry *results) {
    if (root->is_leaf == 1){
        for (int i = 0; i < root->leaf_node.num_entries; i++){
            if (intersects(root->leaf_node.entries[i].mbr, rectangle)){
                results[num_results++] = root->leaf_node.entries[i];
            }
        }
    }
    else{
        for (int i = 0; i < root->non_leaf_node.num_entries; i++){
            if (intersects(root->non_leaf_node.entries[i].mbr, rectangle)){
                searchGetResults(root->non_leaf_node.entries[i].child_ptr, rectangle, results);
            }
        }
    }
}

// Helper Function: Calculates the MBR of a non-leaf node
Rectangle calculateEntryMBR(NonLeafEntry entry){
    // FOR EACH NON LEAF ENTRY; CALCULATE MBR FROM CHILD  NODES
    // FIND -> LOWEST X, LOWEST Y  AND HIGHEST X, HIGHEST Y
    Rectangle mbr;
    NODE next_node = entry.child_ptr;
    int low_x = INFINITY;
    int low_y = INFINITY;
    int high_x;
    int high_y;

    // IF LEAF POINTER; FIND THE COORDINATES FROM ENTRIES
    if (next_node->is_leaf == 1)
    {
        high_x = next_node->leaf_node.entries[0].mbr.top_right.x;
        high_y = next_node->leaf_node.entries[0].mbr.top_right.y;
        for (int i = 0; i < next_node->leaf_node.num_entries; i++)
        {
            Rectangle obj_mbr = next_node->leaf_node.entries[i].mbr;
            low_x = (obj_mbr.bottom_left.x < low_x) ? obj_mbr.bottom_left.x : low_x;
            low_y = (obj_mbr.bottom_left.y < low_y) ? obj_mbr.bottom_left.y : low_y;
            high_x = (obj_mbr.top_right.x > high_x) ? obj_mbr.top_right.x : high_x;
            high_y = (obj_mbr.top_right.y > high_y) ? obj_mbr.top_right.y : high_y;
        }
    }
    else
    {
        // NON LEAF NODE:
        high_x = next_node->non_leaf_node.entries[0].mbr.top_right.x;
        high_y = next_node->non_leaf_node.entries[0].mbr.top_right.y;
        for (int i = 0; i < next_node->non_leaf_node.num_entries; i++)
        {
            Rectangle child_mbr = next_node->non_leaf_node.entries[i].mbr;
            low_x = (child_mbr.bottom_left.x < low_x) ? child_mbr.bottom_left.x : low_x;
            low_y = (child_mbr.bottom_left.y < low_y) ? child_mbr.bottom_left.y : low_y;
            high_x = (child_mbr.top_right.x > high_x) ? child_mbr.top_right.x : high_x;
            high_y = (child_mbr.top_right.y > high_y) ? child_mbr.top_right.y : high_y;
        }
    }
    mbr.bottom_left.x = low_x;
    mbr.bottom_left.y = low_y;
    mbr.top_right.x = high_x;
    mbr.top_right.y = high_y;
    mbr.h = hilbertXYToIndex(5, (mbr.bottom_left.x + mbr.top_right.x) / 2, (mbr.bottom_left.y + mbr.top_right.y) / 2);
    // mbr.h = HilbertValue(mbr.bottom_left, mbr.top_right);
    return mbr;
}

// Helper Funtion: Calculates and Adjusts the LHV values of all non-leaf-nodes in the tree
void adjustLHV(NODE parentNode){
    if (parentNode != NULL)
    {
        for (int i = 0; i < parentNode->non_leaf_node.num_entries; i++)
        {
            parentNode->non_leaf_node.entries[i].largest_hilbert_value = calculateLHV(parentNode->non_leaf_node.entries[i]);
        }
        adjustLHV(parentNode->parent_ptr);
    }
}

// Helper Funtion: Calculates and Adjusts the MBR values of all non-leaf-nodes in the tree
void adjustMBR(NODE parentNode){
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

// Helper Function: Returns a list of siblings of a node with less than MAX_POINTS
//                  Can also be used to find the number of siblings
NODE *cooperatingSiblings(NODE n, int *numSiblings) {
    // taking s=2
    NODE *S = (NODE *)malloc(sizeof(NODE) * MAX_POINTS);

    for (int i = 0; i < MAX_POINTS; i++) {
        S[i] = NULL;
    }

    S[0] = n;
    *numSiblings = 0;
    NODE parentNode = n->parent_ptr;

    if (parentNode == NULL) {
        return S;
    }

    // GO FROM CURRENT NODE TO ROOT FINDING SIBLINGS
    // WITH LESS THAN MAXIMUM POINTERS
    while (parentNode) {
        int index = -1;
        for (int i = 0; i < parentNode->non_leaf_node.num_entries; i++) {
            if (parentNode->non_leaf_node.entries[i].child_ptr == n) {
                index = i;
                break;
            }
        }
        if (index > 0) {
            S[++(*numSiblings)] = parentNode->non_leaf_node.entries[index - 1].child_ptr;
        }
        if (index < parentNode->non_leaf_node.num_entries - 1) {
            S[++(*numSiblings)] = parentNode->non_leaf_node.entries[index + 1].child_ptr;
        }
        n = parentNode;
        parentNode = parentNode->parent_ptr;
    }

    return S;
}


int main(){
    FILE *fp;
    Point points[MAX_POINTS];
    int num_points = 0;

    // Open the file containing the data points
    fp = fopen("data.txt", "r");
    if (fp == NULL)
    {
        printf("Error opening file.\n");
        return 1;
    }
    for (int i = 0; i < MAX_POINTS; i++)
    {
        fscanf(fp, "%d %d\n", &points[num_points].x, &points[num_points].y);
        printf("POINT  = %d %d\n", points[num_points].x, points[num_points].y);
        num_points++;
    }
    fclose(fp);

    // LETS ASSUME THE RECTANGLES INITIALLY ARE THE DATA POINTS
    // AND THEIR CENTRES ARE THE POINTS TOO
    //  CREATE RECTANGLES
    Rectangle rectangles[MAX_POINTS];
    for (int i = 0; i < num_points; i++)
    {
        rectangles[i].bottom_left = points[i];
        rectangles[i].top_right = points[i];
        rectangles[i].h = hilbertXYToIndex(5, points[i].x, points[i].y);
        // rectangles[i].h = HilbertValue(points[i],points[i]);
    }
    // QUICK SORT RECTANGLES ARRAY
    int number_of_rectangles = sizeof(rectangles) / sizeof(rectangles[0]);

    qsort(rectangles, number_of_rectangles, sizeof(struct Rectangle), compare);

    HilbertRTree *Rtree = new_hilbertRTree();
    Rtree->root = new_node(1);
    Rtree->root->parent_ptr = NULL;

    for (int i = 0; i < MAX_POINTS; i++)
    {
        Rtree->root = Insert(Rtree->root, rectangles[i]);
    }
    // traverse_tree(Rtree);
    return 0;
}

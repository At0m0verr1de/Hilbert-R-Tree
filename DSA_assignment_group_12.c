#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_MBR 4
#define MAX_CHILDREN 2
#define BITS_PER_DIM 16

typedef struct Point Point;
struct Point
{
    double x, y;
};

// typedef struct Rectangle Rectangle;
// struct Rectangle
// {
//     Point bottom_left;
//     Point top_right;
//     float h;
// };

// typedef struct LeafEntry LeafEntry;
// struct LeafEntry
// {
//     Rectangle mbr;
//     int obj_id;
// };

// struct LeafNode
// {
//     int num_entries;
//     LeafEntry entries[MAX_CHILDREN];
// };

// typedef struct NonLeafEntry NonLeafEntry;
// struct NonLeafEntry
// {
//     Rectangle mbr;
//     struct Node *child_ptr;
//     int largest_hilbert_value;
// };

// struct NonLeafNode
// {
//     int num_entries;
//     NonLeafEntry entries[MAX_MBR];
// };

// typedef struct Node *NODE;
// struct Node
// {
//     int is_leaf;
//     NODE parent_ptr;
//     union
//     {
//         struct NonLeafNode non_leaf_node;
//         struct LeafNode leaf_node;
//     } u;
// };

// typedef struct HilbertRTree hilbertRTree;
// struct HilbertRTree
// {
//     int height;
//     NODE root;
// };

// /* Sort the entries of a non-leaf node by increasing order of their Hilbert value */
// void sort_entries(NODE n)
// {
//     int i, j;
//     NonLeafNode *nl;
//     LeafNode *ll;
//     NonLeafEntry tmp;
//     LeafEntry ltmp;
//     float hv1, hv2;

//     if (n->is_leaf)
//     {
//         ll = &(n->u.leaf_node);
//         for (i = 0; i < ll->num_entries - 1; i++)
//         {
//             for (j = i + 1; j < ll->num_entries; j++)
//             {
//                 hv1 = hilbert_value(ll->entries[i].mbr.bottom_left);
//                 hv2 = hilbert_value(ll->entries[j].mbr.bottom_left);
//                 if (hv1 > hv2)
//                 {
//                     ltmp = ll->entries[i];
//                     ll->entries[i] = ll->entries[j];
//                     ll->entries[j] = ltmp;
//                 }
//             }
//         }
//     }
//     else
//     {
//         nl = &(n->u.non_leaf_node);
//         for (i = 0; i < nl->num_entries - 1; i++)
//         {
//             for (j = i + 1; j < nl->num_entries; j++)
//             {
//                 hv1 = nl->entries[i].largest_hilbert_value;
//                 hv2 = nl->entries[j].largest_hilbert_value;
//                 if (hv1 > hv2)
//                 {
//                     tmp = nl->entries[i];
//                     nl->entries[i] = nl->entries[j];
//                     nl->entries[j] = tmp;
//                 }
//             }
//         }
//     }
// }

// // Helper function to find the index of the MBR to insert into
// int find_index(NODE n, Rectangle r)
// {
//     int i;
//     for (i = 0; i < n->u.non_leaf_node.num_entries; i++)
//     {
//         if (hilbert_cmp(r.h, n->u.non_leaf_node.entries[i].largest_hilbert_value) < 0)
//             return i;
//     }
//     return n->u.non_leaf_node.num_entries;
// }

// // Helper function to get the MBR of a node
// Rectangle get_mbr(NODE n)
// {
//     Rectangle mbr = n->u.non_leaf_node.entries[0].mbr;
//     int i;
//     for (i = 1; i < n->u.non_leaf_node.num_entries; i++)
//     {
//         mbr = rect_union(mbr, n->u.non_leaf_node.entries[i].mbr);
//     }
//     return mbr;
// }

// // Helper function to adjust the tree after an insertion or deletion
// void adjust_tree(NODE n)
// {
//     if (n == NULL)
//         return;
//     if (n->is_leaf)
//     {
//         // Adjust MBR and LHV in parent levels
//         Rectangle mbr = get_mbr(n);
//         n->parent_ptr->u.non_leaf_node.entries[n->parent_ptr->u.non_leaf_node.num_entries - 1].mbr = mbr;
//         n->parent_ptr->u.non_leaf_node.entries[n->parent_ptr->u.non_leaf_node.num_entries - 1].largest_hilbert_value = n->lhv;
//         adjust_tree(n->parent_ptr);
//     }
//     else
//     {
//         // Adjust MBR and LHV in parent levels
//         Rectangle mbr = get_mbr(n);
//         int i;
//         for (i = 0; i < n->parent_ptr->u.non_leaf_node.num_entries; i++)
//         {
//             if (n->parent_ptr->u.non_leaf_node.entries[i].child_ptr == n)
//             {
//                 n->parent_ptr->u.non_leaf_node.entries[i].mbr = mbr;
//                 n->parent_ptr->u.non_leaf_node.entries[i].largest_hilbert_value = n->lhv;
//                 break;
//             }
//         }
//         adjust_tree(n->parent_ptr);
//     }
// }

// // Helper function to insert an entry into a node
// void insert_entry(NODE n, Rectangle r, int obj_id)
// {
//     if (n->is_leaf)
//     {
//         // Insert the entry into the leaf node
//         n->u.leaf_node.entries[n->u.leaf_node.num_entries].mbr = r;
//         n->u.leaf_node.entries[n->u.leaf_node.num_entries].obj_id = obj_id;
//         n->u.leaf_node.num_entries++;
//         // Adjust the LHV
//         n->lhv = hilbert_value(r.bottom_left);
//         // Adjust the MBR and LHV in parent levels
//         adjust_tree(n);
//     }
//     else
//     {
//         // Find the index of the MBR to insert into
//         int index = find_index(n, r);
//         // Insert the entry into the child node
//         insert_entry(n->u.non_leaf_node.entries[index].child_ptr, r, obj_id);
//     }
// }

// // Helper function to split a leaf node
// void split_leaf(NODE n)
// {
//     int i;
//     // Create a new leaf node
//     NODE new_node = (NODE)malloc(sizeof(struct node));
//     new_node->is_leaf = true;
//     new_node->u.leaf_node.num_entries = 0;
//     new_node->parent_ptr = n->parent_ptr;
//     new_node->lhv = n->lhv;
//     // Copy the last half of the entries to the new node
//     for (i = n->u.leaf_node.num_entries / 2; i < n->u.leaf_node.num_entries; i++)
//     {
//         new_node->u.leaf_node.entries[i - n->u.leaf_node.num_entries / 2] = n->u.leaf_node.entries[i];
//         new_node->u.leaf_node.num_entries++;
//     }
//     // Adjust the number of entries in the old node
//     n->u.leaf_node.num_entries = n->u.leaf_node.num_entries / 2;
//     // Adjust the LHV of the old node
//     n->lhv = hilbert_value(n->u.leaf_node.entries[n->u.leaf_node.num_entries - 1].mbr.bottom_left);
//     // Insert the new node into the parent node
//     insert_entry(n->parent_ptr, new_node->u.leaf_node.entries[0].mbr, new_node->u.leaf_node.entries[0].obj_id);
//     // Adjust the MBR and LHV of the parent node
//     adjust_tree(n->parent_ptr);
// }

// // Helper function to delete an entry from a node
// void delete_entry(NODE n, Rectangle r, int obj_id)
// {
//     // If the node is a leaf node, delete the entry from it
//     if (n->is_leaf)
//     {
//         int i;
//         // Find the index of the entry to delete
//         for (i = 0; i < n->u.leaf_node.num_entries; i++)
//         {
//             if (n->u.leaf_node.entries[i].obj_id == obj_id)
//                 break;
//         }
//         // Shift entries to fill the gap
//         memmove(&n->u.leaf_node.entries[i], &n->u.leaf_node.entries[i + 1], sizeof(LeafEntry) * (n->u.leaf_node.num_entries - i - 1));
//         n->u.leaf_node.num_entries--;
//         // Adjust MBR and LHV in parent levels
//         adjust_tree(n);
//     }
//     else
//     {
//         // If the node is not a leaf node, delete the entry from its child node
//         int index = find_index(n, r);
//         delete_entry(n->u.non_leaf_node.entries[index].child_ptr, r, obj_id);
//     }
// }

// /* Helper function to redistribute entries among siblings */
// void redistribute_entries(NODE n)
// {
//     int i;
//     NODE parent = n->parent_ptr;
//     int index = find_index(parent, get_mbr(n));
//     NODE sibling;
//     int midpoint = parent->u.non_leaf_node.num_entries / 2;
//     int direction;

//     /* Check if n is the leftmost child */
//     if (index == 0)
//     {
//         sibling = parent->u.non_leaf_node.entries[index + 1].child_ptr;
//         direction = 1;
//     }
//     /* Check if n is the rightmost child */
//     else if (index == parent->u.non_leaf_node.num_entries - 1)
//     {
//         sibling = parent->u.non_leaf_node.entries[index - 1].child_ptr;
//         direction = -1;
//     }
//     /* n has two siblings */
//     else
//     {
//         NODE left_sibling = parent->u.non_leaf_node.entries[index - 1].child_ptr;
//         NODE right_sibling = parent->u.non_leaf_node.entries[index + 1].child_ptr;

//         /* Check which sibling has more entries */
//         if (left_sibling->u.non_leaf_node.num_entries > right_sibling->u.non_leaf_node.num_entries)
//         {
//             sibling = left_sibling;
//             direction = -1;
//         }
//         else
//         {
//             sibling = right_sibling;
//             direction = 1;
//         }
//     }

//     /* Check if redistribution is possible */
//     if (sibling->u.non_leaf_node.num_entries > MIN_ENTRIES)
//     {
//         /* Move entries from sibling to n */
//         if (direction == -1)
//         {
//             /* Move rightmost entry from left sibling to n */
//             NODE left_child = sibling->u.non_leaf_node.entries[sibling->u.non_leaf_node.num_entries - 1].child_ptr;
//             Rectangle left_mbr = get_mbr(left_child);
//             parent->u.non_leaf_node.entries[index - 1].mbr = left_mbr;
//             n->u.non_leaf_node.entries[0].mbr = left_mbr;
//             n->u.non_leaf_node.entries[0].child_ptr = left_child;
//             sibling->u.non_leaf_node.num_entries--;
//             n->u.non_leaf_node.num_entries++;
//         }
//         else
//         {
//             /* Move leftmost entry from right sibling to n */
//             NODE right_child = sibling->u.non_leaf_node.entries[0].child_ptr;
//             Rectangle right_mbr = get_mbr(right_child);
//             parent->u.non_leaf_node.entries[index].mbr = right_mbr;
//             n->u.non_leaf_node.entries[n->u.non_leaf_node.num_entries].mbr = right_mbr;
//             n->u.non_leaf_node.entries[n->u.non_leaf_node.num_entries].child_ptr = right_child;
//             sibling->u.non_leaf_node.num_entries--;
//             n->u.non_leaf_node.num_entries++;
//         }
//     }
//     /* Merge with sibling if redistribution is not possible */
//     else
//     {
//         if (direction == -1)
//         {
//             /* Move entries from n to left sibling */
//             for (i = 0; i < n->u.non_leaf_node.num_entries; i++)
//             {
//                 sibling->u.non_leaf_node.entries[sibling->u.non_leaf_node.num_entries + i] = n->u.non_leaf_node.entries[i];
//                 sibling->u.non_leaf_node.num_entries++;
//             }
//             sibling->u.non_leaf_node.entries[sibling->u.non_leaf_node.num_entries].child_ptr = n->u.non_leaf_node.entries[n->u.non_leaf_node.num_entries - 1].child_ptr;
//             sibling->u.non_leaf_node.num_entries++;
//             /* Delete n from parent */
//             delete_entry(parent, get_mbr(n), -1);
//         }
//         else
//         {
//             /* Move entries from n to right sibling */
//             for (i = 0; i < n->u.non_leaf_node.num_entries; i++)
//             {
//                 sibling->u.non_leaf_node.entries[sibling->u.non_leaf_node.num_entries + i + 1] = n->u.non_leaf_node.entries[i];
//                 sibling->u.non_leaf_node.num_entries++;
//             }
//             sibling->u.non_leaf_node.entries[sibling->u.non_leaf_node.num_entries].child_ptr = n->u.non_leaf_node.entries[n->u.non_leaf_node.num_entries - 1].child_ptr;
//             sibling->u.non_leaf_node.num_entries++;
//             /* Delete n from parent */
//             delete_entry(parent, get_mbr(n), -1);
//         }
//     }
// }

// void merge_nodes(NODE n)
// {
//     NODE parent = n->parent_ptr;
//     int index = find_index(parent, get_mbr(n));
//     NODE sibling;
//     if (index == parent->u.non_leaf_node.num_entries - 1)
//     {
//         sibling = parent->u.non_leaf_node.entries[index - 1].child_ptr;
//     }
//     else
//     {
//         sibling = parent->u.non_leaf_node.entries[index + 1].child_ptr;
//     }
//     Rectangle merged_mbr = get_mbr(n);
//     merged_mbr.top_right.x = max(merged_mbr.top_right.x, sibling->u.leaf_node.entries[sibling->u.leaf_node.num_entries - 1].mbr.top_right.x);
//     merged_mbr.top_right.y = max(merged_mbr.top_right.y, sibling->u.leaf_node.entries[sibling->u.leaf_node.num_entries - 1].mbr.top_right.y);
//     merged_mbr.bottom_left.x = min(merged_mbr.bottom_left.x, sibling->u.leaf_node.entries[0].mbr.bottom_left.x);
//     merged_mbr.bottom_left.y = min(merged_mbr.bottom_left.y, sibling->u.leaf_node.entries[0].mbr.bottom_left.y);
//     for (int i = 0; i < sibling->u.leaf_node.num_entries; i++)
//     {
//         insert_entry(n, sibling->u.leaf_node.entries[i].mbr, sibling->u.leaf_node.entries[i].obj_id);
//     }
//     parent->u.non_leaf_node.entries[index].mbr = merged_mbr;
//     parent->u.non_leaf_node.entries[index].largest_hilbert_value = 0;
//     parent->u.non_leaf_node.entries[index].child_ptr = n;
//     parent->u.non_leaf_node.num_entries--;
//     free(sibling);
//     adjust_tree(n);
// }

// Interleave bits of x and y to generate Hilbert index
unsigned int hilbert_index(unsigned int x, unsigned int y)
{
    unsigned int index = 0;
    unsigned int mask = 1 << (BITS_PER_DIM - 1);
    for (int i = 0; i < BITS_PER_DIM; i++)
    {
        unsigned int bit_x = x & mask ? 1 : 0;
        unsigned int bit_y = y & mask ? 1 : 0;
        index |= (bit_x << (2 * i)) | (bit_y << (2 * i + 1));
        mask >>= 1;
    }
    return index;
}

// Calculate Hilbert value of a point
float hilbert_value(Point p)
{
    unsigned int x = (unsigned int)p.x;
    unsigned int y = (unsigned int)p.y;
    unsigned int index = hilbert_index(x, y);
    return (float)index;
}

int main()
{
    Point p;
    p.x = 100;
    p.y = 100;
    printf("%f ", hilbert_value(p));

    return 0;
}
NODE* cooperatingSiblings(NODE n){
    NODE* S = (NODE*)malloc(sizeof(NODE) * MAX_POINTS);
    for (int i = 0; i < MAX_POINTS; i++) {
        S[i] = NULL;
    }
    S[0] = n;
    int numSiblings = 0;
    NODE parentNode = n->parent_ptr;
    //GO FROM CURRENT NODE TO ROOT FINDING SIBLINGS
    //WITH LESS THAN MAXIMUM POINTERS
    while (parentNode) {
            int index = -1;
            for (int i = 0; i < parentNode->u.non_leaf_node.num_entries; i++) {
                if (parentNode->u.non_leaf_node.entries[i].child_ptr == n) {
                    index = i;
                    break;
                }
            }
            if (index > 0 && parentNode->u.non_leaf_node.entries[index - 1].child_ptr->u.leaf_node.num_entries < M) {
                S[++numSiblings] = parentNode->u.non_leaf_node.entries[index - 1].child_ptr;
            } else if (index == 0 && parentNode->u.non_leaf_node.entries[index + 1].child_ptr->u.leaf_node.num_entries < M) {
                S[++numSiblings] = parentNode->u.non_leaf_node.entries[index + 1].child_ptr;
            }
            S[++numSiblings] = parentNode;
            parentNode = parentNode->parent_ptr;
    }
    return S;
}

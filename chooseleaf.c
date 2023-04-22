NODE ChooseLeaf(NODE n, Rectangle r, int h){
    /* RETURNS THE LEAF NODE IN WHICH TO PLACE A NEW RECTANGE*/
    
    /* IF N IS A LEAF, RETURN N*/
    if(n->is_leaf == 1){
        return n;
    }

    /* IF N IS A NON-LEAF NODE, CHOOSE ENTRY (R, PTR, LHV) WITH MINIMUM LHV
    GREATER THAN H*/
    float min_LHV= INFINITY;
    NODE next_node = NULL;
    for(int i = 0; i<n->u.non_leaf_node.num_entries; i++){
            if(n->u.non_leaf_node.entries[i].largest_hilbert_value > h && n->u.non_leaf_node.entries[i].largest_hilbert_value < min_LHV){
                min_LHV = n->u.non_leaf_node.entries[i].largest_hilbert_value;
                next_node = n->u.non_leaf_node.entries[i].child_ptr;
            }
        
    }
    // IF ALL CHILDREN HAVE LHV LESS THAN H
    if(next_node == NULL){
        //CHOOSE THE CHILD NODE WITH LARGEST LHV
        for(int i = 0; i<n->u.non_leaf_node.num_entries; i++){
            
                if(n->u.non_leaf_node.entries[i].largest_hilbert_value > min_LHV){
                    min_LHV = n->u.non_leaf_node.entries[i].largest_hilbert_value;
                    next_node = n->u.non_leaf_node.entries[i].child_ptr;
                }
            
        }
    }

    /* DESCEND UNTIL A LEAF NODE IS REACHED*/
    return ChooseLeaf(next_node, r, h);
}

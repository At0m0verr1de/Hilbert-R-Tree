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
    for(int i = 0; i<MAX_CHILDREN; i++){
        if(n->children[i] != NULL){
            if(n->children[i]->lhv > h && n->children[i]->lhv < min_LHV){
                min_LHV = n->children[i]->lhv;
                next_node = n->children[i];
            }
        }
    }
    // IF ALL CHILDREN HAVE LHV LESS THAN H
    if(next_node == NULL){
        //CHOOSE THE CHILD NODE WITH LARGEST LHV
        for(int i = 0; i<MAX_CHILDREN; i++){
            if(n->children[i] != NULL){
                if(n->children[i]->lhv > min_LHV){
                    min_LHV = n->children[i]->lhv;
                    next_node = n->children[i];
                }
            }
        }
    }

    /* DESCEND UNTIL A LEAF NODE IS REACHED*/
    return ChooseLeaf(next_node, r, h);
}

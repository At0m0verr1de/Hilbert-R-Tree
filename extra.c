NODE createNode(int is_leaf){
    NODE node = (NODE)malloc(sizeof(struct Node));
    
    node->is_leaf = is_leaf;
    node->num_children = 0;

    if(is_leaf == 1){
        node->u.data_points = (int*)malloc(sizeof(int)*M);
        
    }
    else{
        node->u.children = (NODE*)malloc(sizeof(NODE)*M);

        for(int i = 0; i<M; i++){
            node->u.children[i] = NULL;
        }
    }
    return node;

}
void Insert(NODE root, Rectangle r){
    /* INSERTS A NEW RECTANGLE R IN THE HILBERT R-TREE*/
    NODE leafNode = ChooseLeaf(root, r, r.h); //SELECT A LEAF NODE IN WHICH TO PLACE R

    //INSERT R IN A LEAF NODE
        // IF L HAS AN EMPTY SLOT, INSERT R IN L IN APPROPRIATE PLACE ACCORDING TO HILBERT ORDER
        if(leafNode->count < MAX_MBR){
            int i = 0;
            while(i<leafNode->count && leafNode->mbr[i].h < r.h){
                i++;
            }
            //I AT POSITION WHERE MBR[I].H >= R.H
            //SHIFT ALL FORWARD
            for(int k = leafNode->count; j>i; j--){
                leafNode->mbr[j] = leafNode->mbr[j-1];
            }
            //INSERT THE RECTANGLE
            leafNode->mbr[i] = r;
            //INCREMENT THE COUNT OF RECTANGLES IN NODE
            leafNode->count++;
        }
        else{
            // IF L IS FULL: INVOKE HandleOverflow(L, r) RETURNING NEW LEAF
            //IF SPLIT WAS INEVITABLE
            NODE n = HandleOverflow(leafNode, r);
        }
        // PROPOGATE CHANGES UPWARD
        //FORM A SET S CONTAINING L AND ITS COOPERATING SIBLINGS
        //AND NEW LEAF
}

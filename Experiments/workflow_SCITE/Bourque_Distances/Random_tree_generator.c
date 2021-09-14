/*  @copyright: Louxin Zhang, National University of Singapore
 *  This program generates 20000 random labeled  trees of k nodes 
 *  labeled from 1 to n, which are geenrated by applying random
 *  NNI operation on a selected tree in 4000 iterations. At each
 *  iteration, a randome tree is selected from the trees generated
 *  in the previous iteration and NNI operation is appled to it.
 *
 *  In the output file, the trees are listed one by one, in prufer 
 *  tree code format.
 *
 *  Compile command:  gcc Random_tree_generator.c -o generator
 *  Run command: generator (no_of_nodes) <output_tree_file>
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 50 
 
typedef struct {
    unsigned int first;
    unsigned int second;
} edge;
 
edge *prufer_tree(const unsigned int *code, size_t len)
{
    unsigned int i, j, e = 0;
    const unsigned int n = len + 2; /* Nodes */
    edge *tree = malloc((n - 1) * sizeof(edge));
    unsigned int *degrees = calloc(n + 1, sizeof(unsigned int)); /* 1-based */

    if (tree == NULL || degrees == NULL) {
        free(tree);
        free(degrees);
        return NULL;
    }
    /* Set the degrees from the occurrences in the code */
    for (i = 0; i < len; i++) {
        degrees[ code[i]]++;
    }
    /* Add 1 to them all */
    for (i = 1; i <= n; i++) {
        degrees[i]++;
    }
    /* Add edges to nodes in the code */
    for (i = 0; i < len; i++) {
        /* Find the lowest-numbered node with degree 1 */
        for (j = 1; degrees[j] != 1; j++);
        /* Add the edge */
        tree[e].first = j;
        tree[e].second = code[i];
        degrees[j]--;
        degrees[ code[i]]--;
        e++;
    }
    /* Find the last 2 degree-1 nodes */
    for (i = 1; degrees[i] != 1; i++);
    for (j = i + 1; degrees[j] != 1; j++);
    /* Add the last edge */
    tree[e].first = i;
    tree[e].second = j;
    free(degrees);
 
    return tree;
}

int Update(int code[], int size){
 int i;
 int done;

 done=1; 
 for (i=0; i<size-2; i++){
   if (code[i] < size-1) {done=0; break; }
 }
 return done;

}


void Print_Code(FILE *out, int code[], int size){
  int i;

  for (i=0; i<size-2; i++){
    fprintf(out, "%d ", code[i]);
  }
    fprintf(out, "\n");

}

void Add_One_Code(int code[], int size){
   int i;
   int pos;
   int done;

    done =0; pos=size-3; 
    while (done==0){
      if (code[pos]<size-1) { code[pos]=code[pos]+1; done=1;}
      else {
	code[pos]=0; pos=pos-1;
	if (pos==-1) {done=1;}
      }
    }
}

/* from tree to Prufer coder */
void Encode(int tree[][2], int size, int tree_code[]){
 int i, j, k;
 int smalled_leaf;
 int count, s, t;
 int in_node;
 int edge;
 int intree[SIZE][2];
 int degree[SIZE+1];


 for (i=0; i<=size; i++) {
    degree[i]=0;
 }
 for (i=0; i<size-1; i++) {
    intree[i][0]=tree[i][0]; intree[i][1]=tree[i][1];
    degree[tree[i][0]]+=1;
    degree[tree[i][1]]+=1;
 }
 

 s=0;
 t=1;
 while (s<size-2){
  for (i=1; i<=size; i++){
   if (degree[i]==1){
     for (j=0; j<size-1; j++) {
       if (intree[j][0]==i ) { 
	    tree_code[s]=intree[j][1]; 
	    degree[intree[j][1]]=degree[intree[j][1]]-1;
	    intree[j][0]=-1; intree[j][1]=-1;
	    degree[i]=0;  s=1+s;
	    break;
      } 
     if (intree[j][1]==i) {
	    tree_code[s]=intree[j][0]; 
	    degree[intree[j][0]]=degree[intree[j][0]]-1;
	    intree[j][0]=-1; intree[j][1]=-1;
	    degree[i]=0;  s=1+s;
	    break;
     }
    } /* j loop */
    break;
   } 
 } /* i loop */
 } /* while loop */

} /* encode */

/* e=(a, b) so that switch a child of a and b */
void NNInterchange(int in[][2], int out[][2], int size,  int e){
int i, j, k, m;
unsigned int end[2];
unsigned int flag[2];
unsigned int left, right;
unsigned int degree[SIZE+1];

     
     for (j=0; j<=size; j++) { degree[j]=0; }
     for (j=0; j<size-1; j++) {
        out[j][0]=in[j][0]; out[j][1]=in[j][1];
	degree[in[j][0]] +=1; degree[in[j][1]] +=1;
     }
     end[0]=in[e][0]; end[1]=in[e][1];

     if (degree[end[0]]==1) flag[1]=1; else flag[1]=0;
     if (degree[end[1]]==1) flag[0]=1; else flag[0]=0;
     for (m=0; m<size-1; m++) {
         if (m!=e && flag[1]==0){
	   if (out[m][0]==end[0]) {
		out[m][0]=end[1]; flag[1]=1; 
		/* find an edge indicant to end[0] that is not j */
		/* change it to end[1]  */
	   } else if  (out[m][1]==end[0]) {
		out[m][1]=end[1]; flag[1]=1;
		/* find an edge indicant to end[0] that is not j */
		/* change it to end[1]  */
	   }  
         } else if (m!=e && flag[0]==0) {
	   if (out[m][0]==end[1]) {
		out[m][0]=end[0]; flag[0]=1;
	   } else if  (out[m][1]==end[1]) {
		out[m][1]=end[0]; flag[0]=1;
	   }  
         }
	 if ((flag[0]==1) && (flag[1]==1 )) break;
     } /* m loop  */
} /* NNi */


void Print_Tree(edge *tree, int size){
   int i;

   for (i=0; i<size-1; i++) {
     printf("%d %d\n", tree[i].second, tree[i].first);
   }
   printf("\n");
}

void generate_sequence(int size, char *file_name){

unsigned int code[SIZE-2];
unsigned int code2[SIZE-2];
int out_tree[SIZE][2];
int i, j, k, m;
int done;
edge *tree;
/* a tree[i][0] and tree[i][1] are ends of edge  i */
 int int_tree[SIZE][2];
unsigned int tree_code[SIZE];
FILE *out;


done=0;
/*
for (i=0; i<size-2; i++){
  code[i]=0;
  printf("%d ", 0);
}
printf("\n");
*/


 srand(time(0));
  /* generate a random tree */
  for (k=0; k<size-2; k++) { code2[k]=1+rand()%size;  }

  tree=prufer_tree(code2, size-2);  
  Print_Tree(tree, size);

  for (k=0; k<size-1; k++) {
     int_tree[k][0]=tree[k].first; int_tree[k][1]=tree[k].second;
  }
  /* Encode(int_tree, size,  tree_code); Print_Code(tree_code, size); */

  out=fopen(file_name, "w");
  if (out==NULL) { printf("fail to open the output file\n"); exit(100);}

  for (k=0; k<4000;  k++) {
     j=rand()%(size-1);
     NNInterchange(int_tree, out_tree, size,  j);
     Encode(out_tree, size,  tree_code);
     Print_Code(out, tree_code, size);
     /*
     tree=prufer_tree(tree_code, size-2);  Print_Tree(tree, size);
     */

     j=rand()%(size-1);
     NNInterchange(int_tree, out_tree, size, j);
     Encode(out_tree, size, tree_code);
     Print_Code(out, tree_code, size);

     j=rand()%(size-1);
     NNInterchange(int_tree, out_tree, size, j);
     Encode(out_tree, size, tree_code);
     Print_Code(out, tree_code, size);

     j=rand()%(size-1);
     NNInterchange(int_tree, out_tree, size, j);
     Encode(out_tree, size, tree_code);
     Print_Code(out, tree_code, size);

     j=rand()%(size-1);
     NNInterchange(int_tree, out_tree, size,  j);
     Encode(out_tree,size,  tree_code);
     Print_Code(out, tree_code, size);

     for (i=0; i<size-1; i++) {
	int_tree[i][0]=out_tree[i][0]; int_tree[i][1]=out_tree[i][1];
     } 
  } /* k loop */
  fclose(out);

   
}/* generate */


int main(int argc, char *argv[])
{
    unsigned int code[] = {3, 1, 3, 7, 3};
    const size_t len = sizeof(code) / sizeof(unsigned int);
    const unsigned int n = len + 1; /* Edges */
    edge *tree = prufer_tree(code, len);
    unsigned int e;
    int tree0[SIZE][2];
    int tree1[SIZE][2];
    unsigned int tree_code[SIZE];
    int i;
    int size;

    /*
    Encode(tree0, size, tree_code);
    Print_Code(tree_code,  size);

    NNInterchange(tree0,  tree1,  size,  4);
    */

    if (argc !=3) {  
     printf("Command: a.out <no_of_ndoes> <output_tree_file>\n");
     exit(100);
    }
    generate_sequence(atoi(argv[1]), argv[2]);

    return 0;
}


using namespace std;

template<typename MyFloat>
struct grid_struct{

    long grid[3];
    long coords[3];
    MyFloat absval;
    MyFloat delta;

};

template<typename MyFloat>
class Grid{

    //later: think about making some variables/functions private
    //private:
    //int size;

public:
    //constructor:
    Grid(int s=1);
    int size;
    grid_struct<MyFloat>* cells;


    // move constructor
    Grid(Grid &&other) {
        cells = other.cells;
        size = other.size;
        other.cells=NULL;
    }

    ~Grid() {if(cells!=NULL) free(cells);}

    void shift_grid(long s1, long s2, long s3);
    long find_next_ind(long index, int *step);
    long get_index(long x, long y, long z);
    //void check_pbc_grid(long index);
    void check_pbc_grid(long *grid);
    void check_pbc_coords(long index);

};

template<typename MyFloat>
Grid<MyFloat>::Grid(int n){

    cout<< "Constructing grid class" << endl;
    size=n;

    grid_struct<MyFloat> *cells=(grid_struct<MyFloat>*)calloc(n*n*n,sizeof(grid_struct<MyFloat>));
    long g1,g2,g3, gg1,gg2,gg3, ind;
    MyFloat absval;
    for(g1=0;g1<n;g1++){
        gg1=g1;
        if(g1>n/2) gg1=g1-n;
        for(g2=0;g2<n;g2++){
            gg2=g2;
            if(g2>n/2) gg2=g2-n;
            for(g3=0;g3<n;g3++){
                gg3=g3;
                if(g3>n/2) gg3=g3-n;

                ind=(g1*n+g2)*n+g3;
                absval=sqrt(gg1*gg1+gg2*gg2+gg3*gg3);
                cells[ind].absval=absval;
                cells[ind].coords[0]=gg1;
                cells[ind].coords[1]=gg2;
                cells[ind].coords[2]=gg3;
                cells[ind].grid[0]=g1;
                cells[ind].grid[1]=g2;
                cells[ind].grid[2]=g3;
                //if(ind==0){cout<< "in builder "<< cells[0].grid[0] << " "<< cells[0].grid[1] << " "<<cells[0].grid[2]<< endl;} //correct here

            }
        }
    }

    this->cells=cells;

}

template<typename MyFloat>
long Grid<MyFloat>::find_next_ind(long index, int *step){

    long grid[3]={0,0,0};

    grid[0]=this->cells[index].grid[0];//+step[0];
    grid[1]=this->cells[index].grid[1];//+step[1];
    grid[2]=this->cells[index].grid[2];//+step[2];

    // cout << "in finder before sum "<< grid[0] << " " << grid[1] << " "<< grid[2]  <<endl;
    // cout << this->cells[0].grid[0]<< " "<< this->cells[0].grid[1] << " " << this->cells[0].grid[2] <<endl;
    // cout << endl;

    grid[0]=this->cells[index].grid[0]+step[0];
    grid[1]=this->cells[index].grid[1]+step[1];
    grid[2]=this->cells[index].grid[2]+step[2];

    // cout << "in finder after sum "<< grid[0] << " " << grid[1] << " "<< grid[2]  <<endl;
    // cout << this->cells[0].grid[0]<< " "<< this->cells[0].grid[1] << " " << this->cells[0].grid[2] <<endl;
    // cout<< endl;

    this->check_pbc_grid(grid);

    // cout << "after pbc "<< grid[0] << " " << grid[1] << " "<< grid[2]  <<endl;
    // cout << this->cells[0].grid[0]<< " "<< this->cells[0].grid[1] << " " << this->cells[0].grid[2] <<endl;
    // cout<<endl;

    long newind;

    newind=this->get_index(grid[0], grid[1], grid[2]);

    //delete grid;

    return newind;

}

template<typename MyFloat>
void Grid<MyFloat>::check_pbc_grid(long *grid){

    long x=grid[0];
    long y=grid[1];
    long z=grid[2];

    //cout << "in pbc " <<  x << " " << y << " " << z << endl;

    int size=this->size;

    while(x> size-1){x-=size;}
    while(y> size-1){y-=size;}
    while(z> size-1){z-=size;}

    while(x<0){x+=size;}
    while(y<0){y+=size;}
    while(z<0){z+=size;}


    //cout << "in pbc after " <<  x << " " << y << " " << z << endl;

    grid[0]=x;
    grid[1]=y;
    grid[2]=z;

    //cout << "in pbc final" <<  grid[0] << " " << grid[1] << " " << grid[2] << endl;
}

// void Grid::check_pbc_grid_old(long index){

//  long x=this->cells[index].grid[0];
//  long y=this->cells[index].grid[1];
//  long z=this->cells[index].grid[2];

//  int size=this->size;

//   //todo: replace this with single cases with signum or something?

//   while(x> size){x-=size;}
//   while(y> size){y-=size;}
//   while(z> size){z-=size;}

//   while(x<0){x+=size;}
//   while(y<0){y+=size;}
//   while(z<0){z+=size;}

//   this->cells[index].grid[0]=x;
//   this->cells[index].grid[1]=y;
//   this->cells[index].grid[2]=z;

// }

template<typename MyFloat>
void Grid<MyFloat>::check_pbc_coords(long index){

    long x=this->cells[index].coords[0];
    long y=this->cells[index].coords[1];
    long z=this->cells[index].coords[2];

    int size=this->size;

    while(x> size/2){x-=size;}
    while(y> size/2){y-=size;}
    while(z> size/2){z-=size;}

    while(x<-(size/2-1)){x+=size;}
    while(y<-(size/2-1)){y+=size;}
    while(z<-(size/2-1)){z+=size;}

    this->cells[index].coords[0]=x;
    this->cells[index].coords[1]=y;
                this->cells[index].coords[2]=z;


}

template<typename MyFloat>
void Grid<MyFloat>::shift_grid(long s0, long s1, long s2){

    //long coords[3];
    long index;
    long max=(this->size)*(this->size)*(this->size);

    for(index=0; index< max; index++) {
        if(index==0 || index ==1 || index==max-1)
            cout<< "before: index: "<< index << " " << this->cells[index].coords[0] << " " << this->cells[index].coords[1] << " "<< this->cells[index].coords[2] << endl;

        this->cells[index].coords[0]-=s0;
        this->cells[index].coords[1]-=s1;
        this->cells[index].coords[2]-=s2;

        if(index==0 || index ==1 || index==max-1)
            cout<< "intermed: index: "<< index << " " << this->cells[index].coords[0] << " " << this->cells[index].coords[1] << " "<< this->cells[index].coords[2] << endl;

        this->check_pbc_coords(index);

        if(index==0 || index ==1 || index==max-1)
            cout<< "after: index: "<< index << " " << this->cells[index].coords[0] << " " << this->cells[index].coords[1] << " "<< this->cells[index].coords[2] << endl;

    }

    //free(coords);

}

template<typename MyFloat>
long Grid<MyFloat>::get_index(long x, long y, long z){


    long size=this->size;

    // wrap:
    if(x<0) x+=size;
    if(y<0) y+=size;
    if(z<0) z+=size;
    if(x>=size) x-=size;
    if(y>=size) y-=size;
    if(z>=size) z-=size;

    long index=(x*size+y)*size+z;

    return index;

}

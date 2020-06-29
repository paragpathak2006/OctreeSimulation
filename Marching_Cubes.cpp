#include<iostream>
using namespace std;
#include<cstdio>
#include<conio.h>
#define NUMVER 16
#define NUMFACE 90
bool debug=false;
int err=0;
class Face{
public:
	int P[6];
	int n,s;
public:

	Face(){
		for(int j=0;j<6;j++) P[j]=-1;
		n=0;s=0;
	}	
	Face(int i,int p0,int p1,int p2)	{
		for(int j=0;j<6;j++) P[j]=-1;
		if(i==1){	P[0]=p0;		P[1]=p1;		P[2]=p2;}		
		if(i==-1){	P[0]=p2;		P[1]=p1;		P[2]=p0;}		
		n=3;			s=i;
	}	
	Face(int i,int p0,int p1,int p2,int p3)	{
		for(int j=0;j<6;j++) P[j]=-1;
		P[0]=p0;		P[1]=p1;		P[2]=p2;		P[3]=p3;
		n=4;			s=i;
	}	
	Face(int i,int p0,int p1,int p2,int p3,int p4)	{
		for(int j=0;j<6;j++) P[j]=-1;
		P[0]=p0;		P[1]=p1;		P[2]=p2;		P[3]=p3;		P[4]=p4;	
		n=5;			s=i;
	}	
	Face(int i,int p0,int p1,int p2,int p3,int p4,int p5)	{
		for(int j=0;j<6;j++) P[j]=-1;
		P[0]=p0;		P[1]=p1;		P[2]=p2;		P[3]=p3;		P[4]=p4;	P[5]=p5;
		n=6;			s=i;
	}	
	void print(FILE *fp);
};

class Solid{
public:
	Face F[10];
	int n,s;

	int findPos[28];	

	int vertex[NUMVER];
	int numVertex=0;

	int tri[NUMFACE];
	int numTri=0;

public:
	Solid()
	{
		for(int i=0;i<10;i++) F[i]=Face();
		
		for(int i=0;i<NUMFACE;i++) tri[i]=-2;
		for(int i=0;i<NUMVER;i++) vertex[i]=-2;
		
		for(int i=0;i<28;i++) findPos[i]=-5;

		n=0;s=0;	numVertex=0;numTri=0;
	}	
	void add(Face f)	{F[n] = f;		n++;/*for(int i=0;i<6;i++)printf("P[%i]=%i\t",i,f.P[i]);printf("\n");*/}

	void getVertices()	
	{
		for(int i=0;i<10;i++) for(int j=0;j<6;j++) addVertex(F[i].P[j]);
	}	
	void addVertex(int j)	
	{	
		if(j<0) return;
		for(int i=0;i<numVertex;i++) if(vertex[i]==j) return;
		vertex[numVertex]=j;
		findPos[j]=numVertex;
		numVertex++;
		if(numVertex>NUMVER){printf("vertices exeed limits");exit(0);}
	}
	void getTri()	
	{
		for(int i=0;i<10;i++) 
		{
//				printf("numTri=%i\n",numTri);
			if(F[i].n>=3)
			{
				tri[numTri]=findPos[F[i].P[0]];
				tri[numTri+1]=findPos[F[i].P[1]];
				tri[numTri+2]=findPos[F[i].P[2]];
				numTri=numTri+3;
				if(numTri>NUMFACE){printf("Tri exeed limits");exit(0);}
				if(numTri%3!=0){printf("3 numTri is odd %i\n",numTri);getch();}
			}
//				printf("3 numTri=%i\n",numTri);
			if(F[i].n>=4)
			{
				tri[numTri]=findPos[F[i].P[0]];
				tri[numTri+1]=findPos[F[i].P[2]];
				tri[numTri+2]=findPos[F[i].P[3]];
				numTri=numTri+3;
				if(numTri>NUMFACE){printf("Tri exeed limits");exit(0);}
				if(numTri%3!=0){printf("4 numTri is odd %i\n",numTri);getch();}
			}
//				printf("4 numTri=%i\n",numTri);
	
			if(F[i].n>=5)
			{
				tri[numTri]=findPos[F[i].P[0]];
				tri[numTri+1]=findPos[F[i].P[3]];
				tri[numTri+2]=findPos[F[i].P[4]];
				numTri=numTri+3;
				if(numTri>NUMFACE){printf("Tri exeed limits");exit(0);}
				if(numTri%3!=0){printf("5 numTri is odd %i\n",numTri);getch();}
	
			}
//				printf("5 numTri=%i\n",numTri);
	
			if(F[i].n>=6)
			{
				tri[numTri]=findPos[F[i].P[0]];
				tri[numTri+1]=findPos[F[i].P[4]];
				tri[numTri+2]=findPos[F[i].P[5]];
				numTri=numTri+3;
				if(numTri>NUMFACE){printf("Tri exeed limits");exit(0);}
				if(numTri%3!=0){printf("6 numTri is odd %i\n",numTri);getch();}
			}
			//	printf("6 numTri=%i\n",numTri);
		};
	}	
	void setSolid(){getVertices();getTri();}
	void print(FILE *fp);
	void printTri(FILE *fp);	
};

class Cube
{
public:
	int cellType;
	int baseType;
//	int t[13];
	bool set;
	int tP[12];

	Solid pSolid[5];
	int numSolid;
	bool isConjugate;
	bool isAmbiguous;

public:
	void printPoly(FILE *fp);
	void print(FILE *fp);
	void printFV(FILE *fp);
//	void print(FILE *fp,int x,int y);
//	void printBox(FILE *fp,int i,int j);
//	void print(FILE *fp,int sideType,int x,int y);
	
	Cube()
	{
		cellType = -1;
		baseType = -1;
//		for(int i=0;i<13;i++) t[i] = 0;
		for(int i=0;i<12;i++) tP[i] = 0;
		set = false;
		isConjugate = false;
		isAmbiguous = false;
		numSolid=0;
	}

	void setCube(int BaseType,int CellType,bool conjugate){
		isConjugate = conjugate;
		cellType = CellType;
		if(set==true) return;

		baseType = BaseType;
		 for(int i=0;i<5;i++) pSolid[i].setSolid();
		if(BaseType==6||BaseType==3) isAmbiguous=true;
		else isAmbiguous=false;
		set = true;

	}
};

class MarchingCube
{
public:
	Cube cube[256];
	int SPos[28],S[28],P[8],T[13];
	int baseType,cellType,conjugateCellType;

	void setMarchingCube(int type);
//    void getTraingles();
	void getPolyhedras();
public:
	void trimFace(int num,int sense,int i,int j,int k);
	void XP(int num);						//Square
	void XP(int num,char t,int i);			//Traingle or pentagon
	void XP(int num,char t,int i,int j);	//Hexagon
	void XP(int num,char t,char d,char s);	//Rectangle

	void XN(int num);						//Square
	void XN(int num,char t,int i);			//Traingle or pentagon
	void XN(int num,char t,int i,int j);	//Hexagon
	void XN(int num,char t,char d,char s);	//Rectangle

	void YP(int num);						//Square
	void YP(int num,char t,int i);			//Traingle or pentagon
	void YP(int num,char t,int i,int j);	//Hexagon
	void YP(int num,char t,char d,char s);	//Rectangle

	void YN(int num);						//Square
	void YN(int num,char t,int i);			//Traingle or pentagon
	void YN(int num,char t,int i,int j);	//Hexagon
	void YN(int num,char t,char d,char s);	//Rectangle
    
	void ZP(int num);						//Square
	void ZP(int num,char t,int i);			//Traingle or pentagon
	void ZP(int num,char t,int i,int j);	//Hexagon
	void ZP(int num,char t,char d,char s);	//Rectangle

	void ZN(int num);						//Square
	void ZN(int num,char t,int i);			//Traingle or pentagon
	void ZN(int num,char t,int i,int j);	//Hexagon
	void ZN(int num,char t,char d,char s);	//Rectangle


	MarchingCube(){
		
		for(int i=0;i<15;i++)
		{
			setMarchingCube(i);                
			getCellType();
			calculateMarchingCube();
			for(int m=0;m<4;m++)
			{
				for(int l=0;l<4;l++){
			        if(m==0||m==2) Xorient(); 
			        else Zorient(); 
				}
			    Yorient(); 
			}
		
			for(int m=0;m<4;m++)
			{
		    	Zorient();
				for(int l=0;l<4;l++)
			        if(m==0||m==2) Yorient(); 
			}
		}
	}

	void Xorient(){
      forwardFourSwap(P[1],P[5],P[6],P[2]);
        forwardFourSwap(P[0],P[4],P[7],P[3]);

        forwardFourSwap(SPos[21],SPos[25],SPos[26],SPos[22]);
        forwardFourSwap(SPos[20],SPos[24],SPos[27],SPos[23]);
		
        forwardFourSwap(SPos[9],SPos[5],SPos[10],SPos[1]);
        forwardFourSwap(SPos[0],SPos[4],SPos[6],SPos[2]);
        forwardFourSwap(SPos[8],SPos[7],SPos[11],SPos[3]);

        calculateMarchingCube();
	}

	void Yorient(){
        forwardFourSwap(P[7],P[6],P[5],P[4]);
        forwardFourSwap(P[3],P[2],P[1],P[0]);

       forwardFourSwap(SPos[27],SPos[26],SPos[25],SPos[24]);
        forwardFourSwap(SPos[23],SPos[22],SPos[21],SPos[20]);
		
        forwardFourSwap(SPos[7],SPos[6],SPos[5],SPos[4]);
        forwardFourSwap(SPos[11],SPos[10],SPos[9],SPos[8]);
        forwardFourSwap(SPos[3],SPos[2],SPos[1],SPos[0]);

        calculateMarchingCube();
	}

	void Zorient(){
        forwardFourSwap(P[0],P[1],P[5],P[4]);
        forwardFourSwap(P[3],P[2],P[6],P[7]);

        forwardFourSwap(SPos[20],SPos[21],SPos[25],SPos[24]);
        forwardFourSwap(SPos[23],SPos[22],SPos[26],SPos[27]);
		
        forwardFourSwap(SPos[0],SPos[9],SPos[4],SPos[8]);
        forwardFourSwap(SPos[3],SPos[1],SPos[5],SPos[7]);
        forwardFourSwap(SPos[2],SPos[10],SPos[6],SPos[11]);

        calculateMarchingCube();
	}
	void getCellType(){
	int Q[8];
	for(int i=0;i<8;i++) {	if(P[i]==0) Q[i]=1;else Q[i]=0;}	

    cellType =          P[7]*128 + P[6]* 64 + P[5]*32+ P[4]*16 + P[3]*8 + P[2]*4 + P[1]*2 + P[0]*1;   
    conjugateCellType = Q[7]*128 + Q[6]* 64 + Q[5]*32+ Q[4]*16 + Q[3]*8 + Q[2]*4 + Q[1]*2 + Q[0]*1;   

//	    return P[7]*128 + P[6]* 64 + P[5]*32+ P[4]*16 + P[3]*8 + P[2]*4 + P[1]*2 + P[0]*1;   
	} 

	void calculateMarchingCube()
	{
        getCellType();
		for(int i = 0;i<28;i++) S[SPos[i]]=i; 
        if(cube[cellType].set==true) return;

	    getPolyhedras();
		cube[cellType].setCube(baseType,cellType,0);
		cube[conjugateCellType] = cube[cellType];
		cube[conjugateCellType].setCube(baseType,conjugateCellType,1);
	}
    reverseFourSwap(int &i,int &j,int &k,int &l)
	{
        int temp=i;        i=j;        j=k;        k=l;        l=temp;
	}
    forwardFourSwap(int &i,int &j,int &k,int &l)
	{
        int temp=l;        l=k;        k=j;        j=i;        i=temp;
	}
};
int main()
{
	MarchingCube M;
	FILE *fp = fopen("MarchingCubes.txt","w");
//	for(int i = 1;i<101;i++) 	fprintf(fp,"%i\t",i);

	for(int i = 0;i<256;i++) 	M.cube[i].printPoly(fp);
	fclose(fp);

	fp = fopen("Polyhedras.txt","w");
	for(int i = 0;i<256;i++) 	M.cube[i].printFV(fp);
	fclose(fp);

}
	void Cube::print(FILE *fp)
	{
//		fprintf(fp,"%i\t",cellType);		
//		fprintf(fp,"%i\t",baseType);		

//		for(int i=0;i<13;i++) 
//			fprintf(fp,"%i\t",t[i]);

		fprintf(fp,"\n");		
	}
	void Cube::printPoly(FILE *fp)
	{

			if(debug==true)	fprintf(fp,"baseType=%i\tcellType=%i\tisCon=%i\tnumSolids=%i\n",baseType,cellType,isConjugate,numSolid);
		for(int i=0;i<5;i++) 
		{	
			for(int j=0;j<12;j++) fprintf(fp,"%i\t",tP[j]);
//			fprintf(fp,"\t");
//			for(int j=0;j<12;j++) fprintf(fp,"tP(%i)=%i\t",j,tP[j]);

			fprintf(fp,"%i\t%i\t%i\t",isConjugate,isAmbiguous,numSolid);
			pSolid[i].print(fp);
			fprintf(fp,"\n");		
		}
			if(debug==true)	fprintf(fp,"\n");		
	}
	void Cube::printFV(FILE *fp)
	{
		if(debug==true)	fprintf(fp,"baseType=%i\tcellType=%i\tisConju=%i\tnumSolids=%i\n",baseType,cellType,isConjugate,numSolid);
		for(int i=0;i<5;i++) 
		{	
			pSolid[i].printTri(fp);
			fprintf(fp,"\n");		
		}
		if(debug==true)	fprintf(fp,"\n");		

	}

	void Solid::print(FILE *fp)	
	{
		fprintf(fp,"%i\t",n);
		for(int i=0;i<10;i++) 	F[i].print(fp);
	}
	void Solid::printTri(FILE *fp)	
	{
		if(numVertex%2==0) fprintf(fp,"%i\t",numVertex);
		if(numVertex%2==1) fprintf(fp,"**%i**\t",numVertex);
		if(numTri%3==0) fprintf(fp,"%i\t",numTri/3);
		if(numTri%3!=0) fprintf(fp,"**%i**\t",numTri/3);

		if(debug==true)			fprintf(fp,"\t");
		for(int i=0;i<NUMVER;i++) 	fprintf(fp,"%i\t",vertex[i]+1);
		if(debug==true)	fprintf(fp,"\t");
		for(int i=0;i<NUMFACE;i++) 
		{
			if(i%3==0&&debug==true)fprintf(fp,"\t");
			fprintf(fp,"%i\t",tri[i]+1);
		}
		
	}

	void Face::print(FILE *fp)	{
		if(s!=0) 
		{
			fprintf(fp,"%i\t",n);
			for(int i=0;i<6;i++) 
			{
//			printf("\t\t\tPoint(%i)=%i\n",i,P[i]);
			fprintf(fp,"%i\t",P[i]+1);
			}
		}

		if(s==0) 
		{
			fprintf(fp,"%i\t",0);
			for(int i=0;i<6;i++) 
			{
//			printf("\t\t\tPoint(%i)=%i\n",i,P[i]);
			fprintf(fp,"-2\t");
			}
		}
		
		if(s!=-1&&s!=0&&s!=1) {	printf("Print face err s=%i\n",s);getch();}
	}
		void MarchingCube::trimFace(int num,int sense,int i,int j,int k){
			cube[cellType].pSolid[num].add(Face(sense,S[i],S[j],S[k]));		
			cube[cellType].numSolid++;
			if(S[i]>11||S[j]>11||S[k]>11){	printf("err trimFace(%i,%i,%i,%i)\n",S[i],S[j],S[k]);getch();}
			cube[cellType].tP[S[i]]=1;			cube[cellType].tP[S[j]]=1;			cube[cellType].tP[S[k]]=1;
		}
	
		void MarchingCube::XP(int num){
			cube[cellType].pSolid[num].add(Face(1,S[21],S[25],S[26],S[22]));
//				

		}		
		void MarchingCube::XP(int num,char t,int i){
				
			if( t =='T')
			{
				if(i==21) { cube[cellType].pSolid[num].add(Face(1,S[1],S[21],S[9]));return;}
				if(i==25) { cube[cellType].pSolid[num].add(Face(1,S[9],S[25],S[5]));return;}
				if(i==26) { cube[cellType].pSolid[num].add(Face(1,S[5],S[26],S[10]));return;}
				if(i==22) { cube[cellType].pSolid[num].add(Face(1,S[10],S[22],S[1]));return;}
				
			}
			if(t=='P')
			{
				if(i==21) { cube[cellType].pSolid[num].add(Face(1,S[1],S[9],S[25],S[26],S[22]));return;}
				if(i==25) { cube[cellType].pSolid[num].add(Face(1,S[21],S[9],S[5],S[26],S[22]));return;}
				if(i==26) { cube[cellType].pSolid[num].add(Face(1,S[21],S[25],S[5],S[10],S[22]));return;}
				if(i==22) { cube[cellType].pSolid[num].add(Face(1,S[21],S[25],S[26],S[10],S[1]));return;}
			}
			printf("Error in XP(%i,%c,%i)",num,t,i);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
		void MarchingCube::XP(int num,char t,int i,int j){
				

			if(t=='H')
			{
				if(i==21||i==26) { cube[cellType].pSolid[num].add(Face(1,S[1],S[9],S[25],S[5],S[10],S[22]));return;}
				if(i==25||i==22) { cube[cellType].pSolid[num].add(Face(1,S[21],S[9],S[5],S[26],S[10],S[1]));return;}
			}
			printf("Error in XP(%i,%c,%i,%i)",num,t,i,j);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
		void MarchingCube::XP(int num,char t,char d,char s){
				
			if(t=='R')
			{
				if(d=='Y')
				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[10],S[9],S[25],S[26]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[9],S[10],S[22],S[21]));return;}
				}
				if(d=='Z')
				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[1],S[5],S[26],S[22]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[5],S[1],S[21],S[25]));return;}
				}
			}
			printf("Error in XP(%i,%c,%c,%c)",num,t,d,s);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
//---------------------------------------------------------------
		void MarchingCube::XN(int num){
				
			cube[cellType].pSolid[num].add(Face(1,20,23,27,24));
		}		
		void MarchingCube::XN(int num,char t,int i){
				
			if(t=='T')
			{
				if(i==20) { cube[cellType].pSolid[num].add(Face(1,S[8],S[20],S[3]));return;}
				if(i==23) { cube[cellType].pSolid[num].add(Face(1,S[3],S[23],S[11]));return;}
				if(i==27) { cube[cellType].pSolid[num].add(Face(1,S[11],S[27],S[7]));return;}
				if(i==24) { cube[cellType].pSolid[num].add(Face(1,S[7],S[24],S[8]));return;}
			}
			if(t=='P')
			{
				if(i==20) { cube[cellType].pSolid[num].add(Face(1,S[8],S[3],S[23],S[27],S[24]));return;}
				if(i==23) { cube[cellType].pSolid[num].add(Face(1,S[20],S[3],S[11],S[27],S[24]));return;}
				if(i==27) { cube[cellType].pSolid[num].add(Face(1,S[20],S[23],S[11],S[7],S[24]));return;}
				if(i==24) { cube[cellType].pSolid[num].add(Face(1,S[20],S[23],S[27],S[7],S[8]));return;}
			}
			printf("Error in XN(%i,%c,%i)",num,t,i);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
		void MarchingCube::XN(int num,char t,int i,int j)
		{
				
			if(t=='H')
			{
				if(i==20||i==27) { cube[cellType].pSolid[num].add(Face(1,S[8],S[3],S[23],S[11],S[7],S[24]));return;}
				if(i==23||i==24) { cube[cellType].pSolid[num].add(Face(1,S[20],S[3],S[11],S[27],S[7],S[8]));return;}
			}
			printf("Error in XN(%i,%c,%i,%i)",num,t,i,j);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}

		void MarchingCube::XN(int num,char t,char d,char s)
		{
				
			if(t=='R')
			{
				if(d=='Y')
				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[8],S[11],S[27],S[24]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[20],S[23],S[11],S[8]));return;}
				}
				if(d=='Z')
				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[3],S[23],S[27],S[7]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[20],S[3],S[7],S[24]));return;}
				}
			}
			printf("Error in XN(%i,%c,%c,%c)",num,t,d,s);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
//-------------------------------------------------------------------				
	
		void MarchingCube::YP(int num)
		{
				
			cube[cellType].pSolid[num].add(Face(1,S[24],S[27],S[26],S[25]));
		}		
		void MarchingCube::YP(int num,char t,int i)
		{
				
			if(t=='T')
			{
				if(i==24) { cube[cellType].pSolid[num].add(Face(1,S[4],S[24],S[7]));return;}
				if(i==27) { cube[cellType].pSolid[num].add(Face(1,S[7],S[27],S[6]));return;}
				if(i==26) { cube[cellType].pSolid[num].add(Face(1,S[6],S[26],S[5]));return;}
				if(i==25) { cube[cellType].pSolid[num].add(Face(1,S[5],S[25],S[4]));return;}
			}
			if(t=='P')
			{
				if(i==24) { cube[cellType].pSolid[num].add(Face(1,S[4],S[7],S[27],S[26],S[25]));return;}
				if(i==27) { cube[cellType].pSolid[num].add(Face(1,S[24],S[7],S[6],S[26],S[25]));return;}
				if(i==26) { cube[cellType].pSolid[num].add(Face(1,S[24],S[27],S[6],S[5],S[25]));return;}
				if(i==25) { cube[cellType].pSolid[num].add(Face(1,S[24],S[27],S[26],S[5],S[4]));return;}
			}
			printf("Error in YP(%i,%c,%i,%c)",num,t,i);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
		void MarchingCube::YP(int num,char t,int i,int j)
		{
				
			if(t=='H')
			{
				if(i==24||i==26) { cube[cellType].pSolid[num].add(Face(1,S[4],S[7],S[27],S[6],S[5],S[25]));return;}
				if(i==27||i==25) { cube[cellType].pSolid[num].add(Face(1,S[24],S[2],S[6],S[26],S[5],S[4]));return;}
			}
			printf("Error in YP(%i,%c,%i,%i)",num,t,i,j);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
		void MarchingCube::YP(int num,char t,char d,char s)
		{
				
			if(t=='R')
			{
				if(d=='Z')
				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[7],S[27],S[26],S[5]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[24],S[7],S[5],S[25]));return;}
				}
				if(d=='X')
				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[6],S[26],S[25],S[4]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[27],S[6],S[4],S[24]));return;}
				}
			}
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
			printf("Error in YP(%i,%c,%c,%c)",num,t,d,s);	
		}
//---------------------------------------------------------------
		void MarchingCube::YN(int num){
				
			cube[cellType].pSolid[num].add(Face(1,S[20],S[21],S[22],S[23]));
		}		
		void MarchingCube::YN(int num,char t,int i){
				
			if(t=='T')
			{
				if(i==20) { cube[cellType].pSolid[num].add(Face(1,S[3],S[20],S[0]));return;}
				if(i==21) { cube[cellType].pSolid[num].add(Face(1,S[0],S[21],S[1]));return;}
				if(i==22) { cube[cellType].pSolid[num].add(Face(1,S[1],S[22],S[2]));return;}
				if(i==23) { cube[cellType].pSolid[num].add(Face(1,S[2],S[23],S[3]));return;}
			}
			if(t=='P')
			{
				if(i==20) { cube[cellType].pSolid[num].add(Face(1,S[3],S[0],S[21],S[22],S[23]));return;}
				if(i==21) { cube[cellType].pSolid[num].add(Face(1,S[20],S[0],S[1],S[22],S[23]));return;}
				if(i==22) { cube[cellType].pSolid[num].add(Face(1,S[20],S[21],S[1],S[2],S[23]));return;}
				if(i==23) { cube[cellType].pSolid[num].add(Face(1,S[20],S[21],S[22],S[2],S[3]));return;}
			}
			printf("Error in YN(%i,%c,%i)",num,t,i);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
		void MarchingCube::YN(int num,char t,int i,int j){
				
			if(t=='H')
			{
				if(i==20||i==22) { cube[cellType].pSolid[num].add(Face(1,S[3],S[0],S[21],S[1],S[2],S[23]));return;}
				if(i==21||i==23) { cube[cellType].pSolid[num].add(Face(1,S[20],S[0],S[1],S[22],S[2],S[3]));return;}
			}
			printf("Error in YN(%i,%c,%i,%i)",num,t,i,j);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
		void MarchingCube::YN(int num,char t,char d,char s){
				
			if(t=='R')
			{
				if(d=='Z')
				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[1],S[22],S[23],S[3]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[21],S[1],S[3],S[20]));return;}
				}
				if(d=='X')
				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[0],S[21],S[22],S[2]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[20],S[0],S[2],S[23]));return;}
				}
			}
			printf("Error in YN(%i,%c,%c,%c)",num,t,d,s);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
//-------------------------------------------------------------------				
		void MarchingCube::ZP(int num){
				
			cube[cellType].pSolid[num].add(Face(1,S[22],S[26],S[27],S[23]));
		}		
		void MarchingCube::ZP(int num,char t,int i){
				
			if(t=='T')
			{
				if(i==22) { cube[cellType].pSolid[num].add(Face(1,S[2],S[22],S[10]));return;}
				if(i==26) { cube[cellType].pSolid[num].add(Face(1,S[10],S[26],S[6]));return;}
				if(i==27) { cube[cellType].pSolid[num].add(Face(1,S[6],S[27],S[11]));return;}
				if(i==23) { cube[cellType].pSolid[num].add(Face(1,S[11],S[23],S[2]));return;}
			}
			if(t=='P')
			{
				if(i==22) { cube[cellType].pSolid[num].add(Face(1,S[2],S[10],S[26],S[27],S[23]));return;}
				if(i==26) { cube[cellType].pSolid[num].add(Face(1,S[22],S[10],S[6],S[27],S[23]));return;}
				if(i==27) { cube[cellType].pSolid[num].add(Face(1,S[22],S[26],S[6],S[11],S[23]));return;}
				if(i==23) { cube[cellType].pSolid[num].add(Face(1,S[22],S[26],S[27],S[11],S[2]));return;}
			}
				printf("Error in ZP(%i,%c,%i)\n",num,t,i);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
		void MarchingCube::ZP(int num,char t,int i,int j){
				
			if(t=='H')
			{
				if(i==22||i==27) { cube[cellType].pSolid[num].add(Face(1,S[2],S[10],S[26],S[6],S[11],S[23]));return;}
				if(i==26||i==23) { cube[cellType].pSolid[num].add(Face(1,S[22],S[10],S[6],S[27],S[11],S[2]));return;}
				printf("Error in ZP(int num,H,%i,%i)",i,j);	
			}
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
		void MarchingCube::ZP(int num,char t,char d,char s){
				
			if(t=='R')
			{
				if(d=='X')
 				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[2],S[22],S[26],S[6]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[23],S[2],S[6],S[27]));return;}
				}
				if(d=='Y')
				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[10],S[26],S[27],S[11]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[22],S[10],S[11],S[23]));return;}
				}
			}
			printf("Error in ZP(%i,%c,%c,%c)",num,t,d,s);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
//---------------------------------------------------------------
		void MarchingCube::ZN(int num){
				
			cube[cellType].pSolid[num].add(Face(1,S[24],S[25],S[21],S[20]));
		}		
		void MarchingCube::ZN(int num,char t,int i){
				
			if(t=='T')
			{
				if(i==24) { cube[cellType].pSolid[num].add(Face(1,S[8],S[24],S[4]));return;}
				if(i==25) { cube[cellType].pSolid[num].add(Face(1,S[4],S[25],S[9]));return;}
				if(i==21) { cube[cellType].pSolid[num].add(Face(1,S[9],S[21],S[0]));return;}
				if(i==20) { cube[cellType].pSolid[num].add(Face(1,S[0],S[20],S[8]));return;}
			}
			if(t=='P')
			{
				if(i==24) { cube[cellType].pSolid[num].add(Face(1,S[8],S[4],S[25],S[21],S[20]));return;}
				if(i==25) { cube[cellType].pSolid[num].add(Face(1,S[24],S[4],S[9],S[21],S[20]));return;}
				if(i==21) { cube[cellType].pSolid[num].add(Face(1,S[24],S[25],S[9],S[0],S[20]));return;}
				if(i==20) { cube[cellType].pSolid[num].add(Face(1,S[24],S[25],S[21],S[0],S[8]));return;}
			}
				printf("Error in ZN(%i,%c,%i)",num,t,i);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
		void MarchingCube::ZN(int num,char t,int i,int j){
				
			if(t=='H')
			{
				if(i==24||i==21) { cube[cellType].pSolid[num].add(Face(1,S[8],S[4],S[25],S[9],S[0],S[20]));return;}
				if(i==25||i==20) { cube[cellType].pSolid[num].add(Face(1,S[24],S[4],S[9],S[21],S[0],S[8]));return;}
			}
			printf("Error in ZN(%i,%c,%i,%i)\n",num,t,i,j);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
		void MarchingCube::ZN(int num,char t,char d,char s){
				
			if(t=='R')
			{
				if(d=='X')
				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[4],S[25],S[21],S[0]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[24],S[4],S[0],S[20]));return;}
				}
				if(d=='Y')
				{
					if(s=='P') { cube[cellType].pSolid[num].add(Face(1,S[8],S[24],S[25],S[9]));return;}
					if(s=='N') { cube[cellType].pSolid[num].add(Face(1,S[20],S[8],S[9],S[21]));return;}
				}
			}
			printf("Error in ZN(%i,%c,%c,%c)",num,t,d,s);	
			cube[cellType].pSolid[num].add(Face(0,0,0,0));
		}
//-------------------------------------------------------------------				
void MarchingCube::getPolyhedras(){

	if(baseType==0){
//-------------------------------
//-ve Polyhedra
cube[cellType].numSolid = 1;
		XP(0);		XN(0);
		YP(0);		YN(0);
		ZP(0);		ZN(0);
	}

	if(baseType==1){

//-------------------------------
//+ve Polyhedra
		YN(0,'T',20);		XN(0,'T',20);		ZN(0,'T',20);			trimFace(0,1,0,8,3);
//-------------------------------
//-ve Polyhedra
		XP(1);		XN(1,'P',20);
		YP(1);		YN(1,'P',20);
		ZP(1);		ZN(1,'P',20);		
		trimFace(1,-1,0,8,3);
		cube[cellType].numSolid = 2;

	}
	if(baseType==2){
//-------------------------------
//+ve Polyhedra

		XP(0,'T',21);		XN(0,'T',20);
							YN(0,'R','Z','N');
							ZN(0,'R','Y','N');		
		trimFace(0,1,9,8,3);			trimFace(0,1,3,1,9);
//-------------------------------
//-ve Polyhedra
		XP(1,'P',21);		XN(1,'P',20);
		YP(1);				YN(1,'R','Z','P');
		ZP(1);				ZN(1,'R','Y','P');		
		trimFace(1,-1,9,8,3);			trimFace(1,-1,3,1,9);
cube[cellType].numSolid = 2;
	}
	
   	if(baseType==3){
//-------------------------------
//+ve Polyhedra
		YN(0,'T',20);		XN(0,'T',20);		ZN(0,'T',20);			trimFace(0,1,0,8,3);
		XP(1,'T',22);		YN(1,'T',22);		ZP(1,'T',22);			trimFace(1,1,1,2,10);
//-------------------------------
//-ve Polyhedra

		XP(2,'P',22);		XN(2,'P',20);
		YP(2);				YN(2,'H',20,22);
		ZP(2,'P',22);		ZN(2,'P',20);
			trimFace(2,-1,0,8,3);			trimFace(2,-1,1,2,10);
cube[cellType].numSolid = 3;
	}
   	if(baseType==4){

//-------------------------------
//+ve Polyhedra
		XP(0,'R','Z','N');
		YP(0,'R','Z','N');		YN(0,'T',21);
		ZN(0,'P',20);		
		trimFace(0,1,0,1,8);		trimFace(0,1,8,1,7);		trimFace(0,1,7,1,5);
//-------------------------------
//-ve Polyhedra
		XP(1,'R','Z','P');
		YP(1,'R','Z','P');		YN(1,'P',21);
		ZP(1);					ZN(1,'T',20);
		trimFace(1,-1,0,1,8);			trimFace(1,-1,8,1,7);			trimFace(1,-1,7,1,5);
		cube[cellType].numSolid = 2;
	}
   	if(baseType==5){
//-------------------------------
//+ve Polyhedra
	XP(0,'R','Z','N');	XN(0,'R','Z','N');
	YP(0,'R','Z','N');	YN(0,'R','Z','N');
						ZN(0);
	trimFace(0,1,1,5,3);		trimFace(0,1,3,5,7);
//-------------------------------
//-ve Polyhedra

	XP(1,'R','Z','P');	XN(1,'R','Z','P');
	YP(1,'R','Z','P');	YN(1,'R','Z','P');
	ZP(1);
	trimFace(1,-1,1,5,3);		trimFace(1,-1,3,5,7);
cube[cellType].numSolid = 2;
	}
   
   	if(baseType==6){

//-------------------------------
//+ve Polyhedra
	XP(0,'R','Z','N');
	YP(0,'R','Z','N');		YN(0,'T',21);
							ZN(0,'P',20);
	trimFace(0,1,1,5,7);		trimFace(0,1,1,7,8);		trimFace(0,1,1,8,0);
	//-------------------------------
	XN(1,'T',23);		YN(1,'T',23);		ZP(1,'T',23);		trimFace(1,1,1,2,10);

//-------------------------------
//-ve Polyhedra

	XP(2,'R','Z','P');	XN(2,'P',23);
	YP(2,'R','Z','P');	YN(2,'H',21,23);
	ZP(2,'P',23);
	trimFace(2,-1,1,2,10);

	ZN(2,'T',20);		trimFace(2,-1,1,5,7);		trimFace(2,-1,1,7,8);		trimFace(2,-1,1,8,0);
cube[cellType].numSolid = 3;

	}
   	if(baseType==7){

//-------------------------------
//+ve Polyhedra
		XN(0,'T',20);		YN(0,'T',20);		ZN(0,'T',20);			trimFace(0,1,0,8,3);
		XP(1,'T',22);		YN(1,'T',22);		ZP(1,'T',22);			trimFace(1,1,1,2,10);
		XN(2,'T',27);		YP(2,'T',27);		ZP(2,'T',27);			trimFace(2,1,7,6,11);
		XP(3,'T',25);		YP(3,'T',25);		ZN(3,'T',25);			trimFace(3,1,9,5,4);

//-------------------------------
//-ve Polyhedra

		XP(4,'H',22,25);			XN(4,'H',20,27);
		YP(4,'H',25,27);			YN(4,'H',20,22);
		ZP(4,'H',22,27);			ZN(4,'H',20,25);
		trimFace(4,-1,0,8,3);		trimFace(4,-1,1,2,10);			trimFace(4,-1,7,6,11);			trimFace(4,-1,9,5,4);
		cube[cellType].numSolid = 5;
	}
   	if(baseType==8){

//-------------------------------
//+ve Polyhedra

		XP(0,'T',25);		XN(0,'P',23);
		YP(0,'P',26);		YN(0,'T',20);	
		ZP(0,'T',27);		ZN(0,'P',21);
		trimFace(0,1,5,0,9);		trimFace(0,1,6,0,5);		trimFace(0,1,6,3,0);		trimFace(0,1,11,3,6);

//-------------------------------
//-ve Polyhedra

		XP(1,'P',25);		XN(1,'T',23);
		YP(1,'T',26);		YN(1,'P',20);	
		ZP(1,'P',27);		ZN(1,'T',21);
		trimFace(1,-1,5,0,9);		trimFace(1,-1,6,0,5);		trimFace(1,-1,6,3,0);		trimFace(1,-1,11,3,6);
cube[cellType].numSolid = 2;
}
   	if(baseType==9){

//-------------------------------
//+ve Polyhedra
	XP(0,'R','Z','N');			XN(0,'R','Y','P');
	YP(0,'P',26);				YN(0,'T',21);
	ZP(0,'T',27);				ZN(0,'P',20);
	trimFace(0,1,1,5,0);		trimFace(0,1,0,5,11);		trimFace(0,1,11,5,6);			trimFace(0,1,0,11,8);
	
//-------------------------------
//-ve Polyhedra
	
	XP(1,'R','Z','P');	XN(1,'R','Y','N');
	YP(1,'T',26);			YN(1,'P',21);
	ZP(1,'P',27);			ZN(1,'T',20);
	trimFace(1,-1,1,5,0);		trimFace(1,-1,0,5,11);	trimFace(1,-1,11,5,6);		trimFace(1,-1,0,11,8);
cube[cellType].numSolid = 2;
}
   	if(baseType==10){

//-------------------------------
//+ve Polyhedra
		YN(0,'T',20);		XN(0,'T',20);		ZN(0,'T',20);			trimFace(0,1,0,8,3);
		YP(0,'T',26);		XP(0,'T',26);		ZP(0,'T',26);			trimFace(0,1,6,5,10);
//-------------------------------
//-ve Polyhedra
		XP(1,'P',26);		XN(1,'P',20);
		YP(1,'P',26);		YN(1,'P',20);
		ZP(1,'P',26);		ZN(1,'P',20);
		trimFace(1,-1,6,5,10);			trimFace(1,-1,0,8,3);
		cube[cellType].numSolid = 2;
	}
   	if(baseType==11){

//-------------------------------
//+ve Polyhedra

		XP(0,'T',21);		XN(0,'T',20);
						YN(0,'R','Z','N');
						ZN(0,'R','Y','N');		
		trimFace(0,1,9,8,3);			trimFace(0,1,3,1,9);

		YP(1,'T',26);		XP(1,'T',26);		ZP(1,'T',26);			trimFace(1,1,6,5,10);
//-------------------------------
//-ve Polyhedra

		XP(2,'H',21,26);		XN(2,'P',20);
		YP(2,'P',26);			YN(2,'R','Z','P');	
		ZP(2,'P',26);			ZN(2,'R','Y','P');
		trimFace(2,-1,9,8,3);			trimFace(2,-1,3,1,9);
cube[cellType].numSolid = 3;
	}
   	if(baseType==12){
//-------------------------------
//+ve Polyhedra

		YP(0,'T',26);		XP(0,'T',26);		ZP(0,'T',26);			trimFace(0,1,6,5,10);
		YN(1,'T',23);		XN(1,'T',23);		ZP(1,'T',23);			trimFace(1,1,2,3,11);
		YN(2,'T',21);		XP(2,'T',21);		ZN(2,'T',21);			trimFace(2,1,0,1,9);
//-------------------------------
//-ve Polyhedra
		XP(3,'H',26,21);		XN(3,'P',23);
		YP(3,'P',26);			YN(3,'H',23,21);
		ZP(3,'H',26,23);		ZN(3,'P',21);
		trimFace(3,-1,2,3,11);			trimFace(3,-1,6,5,10);			trimFace(3,-1,0,1,9);

		cube[cellType].numSolid = 4;
	}
   	if(baseType==13){
//-------------------------------
		//+ve Polyhedra
		XP(0,'R','Y','P');
		YP(0,'R','X','P');
		ZP(0,'T',26);			ZN(0,'T',25);
		trimFace(0,1,6,9,10);		trimFace(0,1,6,4,9);
		//-------------------------------
		
						XN(1,'R','Y','N');
						YN(1,'R','X','N');
		ZP(1,'T',23);		ZN(1,'T',20);
		trimFace(1,1,0,8,11);		trimFace(1,1,0,11,2);
//-------------------------------
//-ve Polyhedra

		XP(2,'R','Y','N');		XN(2,'R','Y','P');
		YP(2,'R','X','N');		YN(2,'R','X','P');
		ZP(2,'H',26,23);			ZN(2,'H',25,20);
		trimFace(2,-1,0,8,11);			trimFace(2,-1,0,11,2);			trimFace(2,-1,6,9,10);			trimFace(2,-1,6,4,9);

		cube[cellType].numSolid = 3;
	}
   	if(baseType==14){
		
//-------------------------------
//+ve Polyhedra
		XP(0,'R','Y','P');	XN(0,'R','Z','N');
		YP(0,'P',27);			YN(0,'T',20);
		ZP(0,'T',26);			ZN(0,'P',21);
		trimFace(0,1,0,9,10);		trimFace(0,1,7,10,6);		trimFace(0,1,3,0,7);		trimFace(0,1,0,10,7);

//-------------------------------
//-ve Polyhedra
		XP(1,'R','Y','N');	XN(1,'R','Z','P');
		YP(1,'T',27);			YN(1,'P',20);
		ZP(1,'P',26);			ZN(1,'T',21);
		trimFace(1,-1,0,9,10);		trimFace(1,-1,7,10,6);		trimFace(1,-1,3,0,7);		trimFace(1,-1,0,10,7);

		cube[cellType].numSolid = 2;
	}
}


void MarchingCube::setMarchingCube(int type)
{
	baseType = type;
	cellType = -1;

	for(int i=0;i<8;i++) P[i]=0;		
	for(int i=0;i<28;i++) SPos[i]=i;		
	for(int i=0;i<28;i++) S[i]=i;		
	for(int i=0;i<13;i++) T[i]=0;

	if(type==1){
		P[0]=1;
	}
	if(type==2){
		P[0]=1;		P[1]=1;
	}
	if(type==3){
		P[0]=1;		P[2]=1;
	}
	if(type==4)	{
		P[1]=1;		P[5]=1;		P[4]=1;
	}
	if(type==5) {
		P[1]=1;		P[0]=1;		P[4]=1;		P[5]=1;
	}
	if(type==6)	{
		P[1]=1;		P[3]=1;		P[4]=1;		P[5]=1;
	}
	if(type==7)	{
		P[2]=1;		P[0]=1;		P[7]=1;		P[5]=1;
	}
	if(type==8)	{
		P[7]=1;		P[0]=1;		P[4]=1;		P[5]=1;
	}
	if(type==9)	{
		P[1]=1;		P[7]=1;		P[4]=1;		P[5]=1;
	}
	if(type==10)	{
		P[0]=1;		P[6]=1;
	}
	if(type==11)	{
		P[1]=1;		P[0]=1;		P[6]=1;
	}
	if(type==12)	{
		P[1]=1;		P[3]=1;		P[6]=1;
	}
	if(type==13)	{
		P[3]=1;		P[0]=1;		P[6]=1;		P[5]=1;
	}
	if(type==14)	{
		P[6]=1;		P[0]=1;		P[4]=1;		P[5]=1;
	}
}             

//void MarchingCube::getTraingles(){
//
//	if(baseType==0){
////-------------------------------
////-ve Polyhedra
//	}
//
//	if(baseType==1){
//		T[0]= 1;
//		T[1]=S[0];
//		T[2]=S[3];
//		T[3]=S[8];
//
//
//
//	}
//	if(baseType==2){
//		T[0]= 2;
//		T[1]=S[8];    		T[2]=S[9];    		T[3]=S[3];
//		T[4]=S[3];    		T[5]=S[1];    		T[6]=S[9];
//
//	}
//	
//   	if(baseType==3){
//		T[0]= 2;
//		T[1]=S[3];    		T[2]=S[8];    		T[3]=S[0];
//		T[4]=S[2];    		T[5]=S[1];    		T[6]=S[10];
//	}
//   	if(baseType==4){
//		T[0]= 3;
//		T[1]=S[8];    		T[2]=S[0];    		T[3]=S[1];
//		T[4]=S[1];    		T[5]=S[7];    		T[6]=S[8];
//		T[7]=S[1];    		T[8]=S[7];    		T[9]=S[5];
//
//	}
//   	if(baseType==5){
//		T[0]= 2;
//		T[1]=S[3];    		T[2]=S[7];    		T[3]=S[5];
//		T[4]=S[3];    		T[5]=S[1];    		T[6]=S[5];
//
//	}
//   	if(baseType==6){
//		T[0]= 3;
//		T[1]=S[3];    		T[2]=S[2];    		T[3]=S[11];
//		T[4]=S[1];    		T[5]=S[7];    		T[6]=S[8];
//		T[7]=S[1];    		T[8]=S[7];    		T[9]=S[5];
//
//	}
//   	if(baseType==7){
//		T[0]= 4;
//		T[1]=S[3];    		T[2]=S[8];    		T[3]=S[0];
//		T[4]=S[4];    		T[5]=S[5];    		T[6]=S[9];
//		T[7]=S[1];    		T[8]=S[2];    		T[9]=S[10];
//		T[10]=S[11];   		T[11]=S[7];    		T[12]=S[6];
//	}
//   	if(baseType==8){
//		T[0]= 4;
//		T[1]=S[3];    		T[2]=S[11];    		T[3]=S[6];
//		T[4]=S[3];    		T[5]=S[0];    		T[6]=S[6];
//		T[7]=S[5];    		T[8]=S[6];    		T[9]=S[0];
//		T[10]=S[5];    		T[11]=S[9];    		T[12]=S[0];
//
//}
//   	if(baseType==9){
//		T[0]= 4;
//		T[1]=S[11];    		T[2]=S[8];    		T[3]=S[0];
//		T[4]=S[0];    		T[5]=S[1];    		T[6]=S[5];
//		T[7]=S[5];    		T[8]=S[6];    		T[9]=S[11];
//		T[10]=S[11];    	T[11]=S[5];    		T[12]=S[0];
//
//}
//   	if(baseType==10){
//		T[0]= 2;
//		T[1]=S[3];    		T[2]=S[8];    		T[3]=S[0];
//		T[4]=S[10];    		T[5]=S[6];    		T[6]=S[5];
//
//	}
//   	if(baseType==11){
//		T[0]= 3;
//		T[1]=S[3];    		T[2]=S[8];    		T[3]=S[9];
//		T[4]=S[1];    		T[5]=S[9];    		T[6]=S[3];
//		T[7]=S[10];    		T[8]=S[6];    		T[9]=S[5];
//
//	}
//   	if(baseType==12){
//		T[0]= 3;
//		T[1]=S[3];    		T[2]=S[7];    		T[3]=S[11];
//		T[4]=S[1];    		T[5]=S[0];    		T[6]=S[9];
//		T[7]=S[6];    		T[8]=S[10];    		T[9]=S[5];
//	}
//   	if(baseType==13){
//		T[0]= 4;
//		T[1]=S[11];    		T[2]=S[8];    		T[3]=S[0];
//		T[4]=S[11];    		T[5]=S[8];    		T[6]=S[2];
//		T[7]=S[4];    		T[8]=S[6];    		T[9]=S[9];
//		T[10]=S[10];    	T[11]=S[6];    		T[12]=S[9];
//	}
//   	if(baseType==14){
//		T[0]= 4;
//		T[1]=S[3];    		T[2]=S[7];    		T[3]=S[0];
//		T[4]=S[0];    		T[5]=S[9];    		T[6]=S[10];
//		T[7]=S[10];    		T[8]=S[6];    		T[9]=S[7];
//		T[10]=S[10];    	T[11]=S[7];    		T[12]=S[0];
//		
//	}
//}
//
//
//	void Cube::print(FILE *fp,int x,int y)
//	{
//		printBox(fp,x,y);
//		for(int i=0;i<t[0];i++) {
//			print(fp,t[3*i+1],x,y);
//			print(fp,t[3*i+2],x,y);
//			print(fp,t[3*i+3],x,y);
//			print(fp,t[3*i+1],x,y);
//		}
//		fprintf(fp,"\n\n\n\n");		
//	}
//
//	void Cube::printBox(FILE *fp,int i,int j)
//	{
//		int x=3*i;
//		int y=3*j;
//		int z=0;//3*j;
//		fprintf(fp,"%i\t%i\t%i\n",-1+x,-1+y,-1+z);
//		fprintf(fp,"%i\t%i\t%i\n",1+x,-1+y,-1+z);
//		fprintf(fp,"%i\t%i\t%i\n",1+x,-1+y,1+z);
//		fprintf(fp,"%i\t%i\t%i\n",-1+x,-1+y,1+z);
//		fprintf(fp,"%i\t%i\t%i\n",-1+x,-1+y,-1+z);
//
////		fprintf(fp,"\n\n\n\n");		
//
//		fprintf(fp,"%i\t%i\t%i\n",-1+x,1+y,-1+z);
//		fprintf(fp,"%i\t%i\t%i\n",1+x,1+y,-1+z);
//		fprintf(fp,"%i\t%i\t%i\n",1+x,1+y,1+z);
//		fprintf(fp,"%i\t%i\t%i\n",-1+x,1+y,1+z);
//		fprintf(fp,"%i\t%i\t%i\n",-1+x,1+y,-1+z);
//		
//		fprintf(fp,"\n\n\n\n");		
//
//		fprintf(fp,"%i\t%i\t%i\n",1+x,-1+y,-1+z);
//		fprintf(fp,"%i\t%i\t%i\n",1+x,1+y,-1+z);
//
//		fprintf(fp,"\n\n\n\n");		
//
//		fprintf(fp,"%i\t%i\t%i\n",1+x,-1+y,1+z);
//		fprintf(fp,"%i\t%i\t%i\n",1+x,1+y,1+z);
//
//		fprintf(fp,"\n\n\n\n");		
//
//		fprintf(fp,"%i\t%i\t%i\n",-1+x,-1+y,1+z);
//		fprintf(fp,"%i\t%i\t%i\n",-1+x,1+y,1+z);
//
//		fprintf(fp,"\n\n\n\n");		
//
//
//	}
//
//	void Cube::print(FILE *fp,int sideType,int x,int y)
//	{
//		int i=3*x;
//		int j=3*y;
//		int k =0;// 3*y;
//		
//		switch(sideType){
//			case 0:	fprintf(fp,"%i\t%i\t%i\n",0+i,-1+j,-1+k);			break;
//			case 1:	fprintf(fp,"%i\t%i\t%i\n",1+i,-1+j,0+k);			break;
//			case 2:	fprintf(fp,"%i\t%i\t%i\n",0+i,-1+j,1+k);			break;
//			case 3:	fprintf(fp,"%i\t%i\t%i\n",-1+i,-1+j,0+k);			break;
//
//			case 4:	fprintf(fp,"%i\t%i\t%i\n",0+i,1+j,-1+k);			break;
//			case 5:	fprintf(fp,"%i\t%i\t%i\n",1+i,1+j,0+k);			break;
//			case 6:	fprintf(fp,"%i\t%i\t%i\n",0+i,1+j,1+k);			break;
//			case 7:	fprintf(fp,"%i\t%i\t%i\n",-1+i,1+j,0+k);			break;
//
//			case 8:	fprintf(fp,"%i\t%i\t%i\n",-1+i,0+j,-1+k);			break;
//			case 9:	fprintf(fp,"%i\t%i\t%i\n",1+i,0+j,-1+k);			break;
//			case 10:fprintf(fp,"%i\t%i\t%i\n",1+i,0+j,1+k);			break;
//			case 11:fprintf(fp,"%i\t%i\t%i\n",-1+i,0+j,1+k);			break;
//		}		
//	}


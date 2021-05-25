#ifndef _FMM2DTree_HPP__
#define _FMM2DTree_HPP__

#include <vector>
#include <Eigen/Dense>

#define EIGEN_DONT_PARALLELIZE

const double PI	=	3.1415926535897932384;

struct pts2D {
	double x,y;
};

class FMM2DBox {
public:
	int boxNumber;
	int parentNumber;
	int childrenNumbers[4];
	int neighborNumbers[8];
	int innerNumbers[16];
	int outerNumbers[24];

	FMM2DBox () {
		boxNumber		=	-1;
		parentNumber	=	-1;
		for (int l=0; l<4; ++l) {
			childrenNumbers[l]	=	-1;
		}
		for (int l=0; l<8; ++l) {
			neighborNumbers[l]	=	-1;
		}
		for (int l=0; l<16; ++l) {
			innerNumbers[l]		=	-1;
		}
		for (int l=0; l<24; ++l) {
			outerNumbers[l]		=	-1;
		}
	}

	Eigen::VectorXd multipoles;
	Eigen::VectorXd locals;

	pts2D center;

	//	The following will be stored only at the leaf nodes
	std::vector<pts2D> chebNodes;
};

class kernel {
public:
	bool isTrans;		//	Checks if the kernel is translation invariant, i.e., the kernel is K(r).
	bool isHomog;		//	Checks if the kernel is homogeneous, i.e., K(r) = r^{alpha}.
	bool isLogHomog;	//	Checks if the kernel is log-homogeneous, i.e., K(r) = log(r^{alpha}).
	double alpha;		//	Degree of homogeneity of the kernel.
	kernel() {};
	~kernel() {};
	virtual double getInteraction(const pts2D r1, const pts2D r2, double a){
		return 0.0;
	};	//	Kernel entry generator
};

template <typename kerneltype>
class FMM2DTree {
public:
	kerneltype* K;
	int nLevels;			//	Number of levels in the tree.
	int nChebNodes;			//	Number of Chebyshev nodes along one direction.
	int rank;				//	Rank of interaction, i.e., rank = nChebNodes*nChebNodes.
	int N;					//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	double smallestBoxSize;	//	This is L/2.0^(nLevels).
	double a;				//	Cut-off for self-interaction. This is less than the length of the smallest box size.

	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> boxRadius;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<double> boxHomogRadius;			//	Stores the value of boxRadius^{alpha}
	std::vector<double> boxLogHomogRadius;		//	Stores the value of alpha*log(boxRadius)
	std::vector<std::vector<FMM2DBox> > tree;	//	The tree storing all the information.

	//	Chebyshev nodes
	std::vector<double> standardChebNodes1D;
	std::vector<pts2D> standardChebNodes;
	std::vector<pts2D> standardChebNodesChild;
	std::vector<pts2D> leafChebNodes;

	//	Different Operators
	Eigen::MatrixXd selfInteraction;		//	Needed only at the leaf level.
	Eigen::MatrixXd neighborInteraction[8];	//	Neighbor interaction only needed at the leaf level.
	Eigen::MatrixXd M2M[4];					//	Transfer from multipoles of 4 children to multipoles of parent.
	Eigen::MatrixXd L2L[4];					//	Transfer from locals of parent to locals of 4 children.
	Eigen::MatrixXd M2LInner[16];			//	M2L of inner interactions. This is done on the box [-L,L]^2.
	Eigen::MatrixXd M2LOuter[24];			//	M2L of outer interactions. This is done on the box [-L,L]^2.

// public:
	FMM2DTree(kerneltype* K, int nLevels, int nChebNodes, double L) {
		this->K					=	K;
		this->nLevels			=	nLevels;
		this->nChebNodes		=	nChebNodes;
		this->rank				=	nChebNodes*nChebNodes;
		this->L					=	L;
		nBoxesPerLevel.push_back(1);
		boxRadius.push_back(L);
		boxHomogRadius.push_back(pow(L,K->alpha));
		boxLogHomogRadius.push_back(K->alpha*log(L));
		for (int k=1; k<=nLevels; ++k) {
			nBoxesPerLevel.push_back(4*nBoxesPerLevel[k-1]);
			boxRadius.push_back(0.5*boxRadius[k-1]);
			boxHomogRadius.push_back(pow(0.5,K->alpha)*boxHomogRadius[k-1]);
			boxLogHomogRadius.push_back(boxLogHomogRadius[k-1]-K->alpha*log(2));
		}
		this->smallestBoxSize	=	boxRadius[nLevels];
		this->a					=	smallestBoxSize;
		this->N					=	rank*nBoxesPerLevel[nLevels];
	}

	std::vector<pts2D> shift_Cheb_Nodes(double xShift, double yShift) {
		std::vector<pts2D> shiftedChebNodes;
		for (int k=0; k<rank; ++k) {
			pts2D temp;
			temp.x	=	standardChebNodes[k].x+2*xShift;
			temp.y	=	standardChebNodes[k].y+2*yShift;
			shiftedChebNodes.push_back(temp);
		}
		return shiftedChebNodes;
	}
	
	std::vector<pts2D> shift_Leaf_Cheb_Nodes(double xShift, double yShift) {
		std::vector<pts2D> shiftedChebNodes;
		for (int k=0; k<rank; ++k) {
			pts2D temp;
			temp.x	=	leafChebNodes[k].x+xShift;
			temp.y	=	leafChebNodes[k].y+yShift;
			shiftedChebNodes.push_back(temp);
		}
		return shiftedChebNodes;
	}

	//	get_ChebPoly
	double get_ChebPoly(double x, int n) {
		return cos(n*acos(x));
	}

	//	get_S
	double get_S(double x, double y, int n) {
		double S	=	0.5;
		for (int k=1; k<n; ++k) {
			S+=get_ChebPoly(x,k)*get_ChebPoly(y,k);
		}
		return 2.0/n*S;
	}
	//	set_Standard_Cheb_Nodes
	void set_Standard_Cheb_Nodes() {
		for (int k=0; k<nChebNodes; ++k) {
			standardChebNodes1D.push_back(-cos((k+0.5)/nChebNodes*PI));
		}
		pts2D temp1;
		for (int j=0; j<nChebNodes; ++j) {
			for (int k=0; k<nChebNodes; ++k) {
				temp1.x	=	standardChebNodes1D[k];
				temp1.y	=	standardChebNodes1D[j];
				standardChebNodes.push_back(temp1);
			}
		}
		//	Left Bottom child, i.e., Child 0
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x-0.5;
				temp1.y	=	0.5*temp1.y-0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Right Bottom child, i.e., Child 1
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x+0.5;
				temp1.y	=	0.5*temp1.y-0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Right Top child, i.e., Child 2
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x+0.5;
				temp1.y	=	0.5*temp1.y+0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Left Top child, i.e., Child 3
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x-0.5;
				temp1.y	=	0.5*temp1.y+0.5;
				standardChebNodesChild.push_back(temp1);
		}
	}

	void get_Transfer_Matrix() {
		for (int l=0; l<4; ++l) {
			L2L[l]	=	Eigen::MatrixXd(rank,rank);
			for (int j=0; j<rank; ++j) {
				for (int k=0; k<rank; ++k) {
					L2L[l](j,k)	=	get_S(standardChebNodes[k].x, standardChebNodesChild[j+l*rank].x, nChebNodes)*get_S(standardChebNodes[k].y, standardChebNodesChild[j+l*rank].y, nChebNodes);
				}
			}
		}
		for (int l=0; l<4; ++l) {
			M2M[l]	=	L2L[l].transpose();
		}
	}

	void createTree() {
		//	First create root and add to tree
		FMM2DBox root;
		root.boxNumber		=	0;
		root.parentNumber	=	-1;
		#pragma omp parallel for
		for (int l=0; l<4; ++l) {
			root.childrenNumbers[l]	=	l;
		}
		#pragma omp parallel for
		for (int l=0; l<8; ++l) {
			root.neighborNumbers[l]	=	-1;
		}
		#pragma omp parallel for
		for (int l=0; l<16; ++l) {
			root.innerNumbers[l]	=	-1;
		}
		#pragma omp parallel for
		for (int l=0; l<24; ++l) {
			root.outerNumbers[l]	=	-1;
		}
		std::vector<FMM2DBox> rootLevel;
		rootLevel.push_back(root);
		tree.push_back(rootLevel);

		for (int j=1; j<=nLevels; ++j) {
			std::vector<FMM2DBox> level;
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				FMM2DBox box;
				box.boxNumber		=	k;
				box.parentNumber	=	k/4;
				for (int l=0; l<4; ++l) {
					box.childrenNumbers[l]	=	4*k+l;
				}
				level.push_back(box);
			}
			tree.push_back(level);
		}
	}

	//	Assigns the interactions for child0 of a box
	void assign_Child0_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N5  |  N4  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  **  |  N3  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[3]	=	nC+1;
			tree[nL][nC].neighborNumbers[4]	=	nC+2;
			tree[nL][nC].neighborNumbers[5]	=	nC+3;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/****************************/
			/*				   ______	*/
			/*				  |		 |	*/
			/*				  |	 **  |	*/
			/*	 _____________|______|  */
			/*	|	   |	  |			*/
			/*	|  I15 |  N0  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I0  |  I1  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 ______			  		*/
			/*	|	   |	  			*/
			/*	|  **  |				*/
			/*	|______|______			*/
			/*	|	   |	  |			*/
			/*	|  N1  |  N2  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I2  |  I3  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's second neighbor
		{
			/************************************/
			/*	 ______			  				*/
			/*	|	   |	  					*/
			/*	|  **  |						*/
			/*	|______|	   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  I5  |  O8  |	*/
			/*				  |______|______|	*/
			/*				  |	     |	    |	*/
			/*				  |  I4  |  O7  |	*/
			/*				  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[8]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's third neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  I7  |  O10 |	*/
			/*	 ______		  |______|______|	*/
			/*	|	   |	  |	     |	    |	*/
			/*	|  **  |	  |  I6  |  O9  |	*/
			/*	|______|	  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN!=-1) {
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[9]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[3];
			}

		}

		//	Assign children of parent's fourth neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  O13 |  O12 |	*/
			/*				  |______|______|	*/
			/*				  |	     |	    |	*/
			/*		    	  |  I8  |  O11 |	*/
			/*		    	  |______|______|	*/
			/*									*/
			/*									*/
			/*	 ______							*/
			/*  |      |						*/
			/*  |  **  |						*/
			/*  |______|						*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[11]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[12]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[13]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  O15 |  O14 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  I10 |  I9  |		*/
			/*	|______|______|		*/
			/*						*/
			/*						*/
			/*	 ______				*/
			/*  |	   |			*/
			/*	|  **  |			*/
			/*	|______|			*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[14]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[15]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/****************************/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  O17 |  O16 |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I12 |  I11 |			*/
			/*	|______|______|			*/
			/*							*/
			/*							*/
			/*				   ______	*/
			/*  			  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[16]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[17]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/****************************/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I13 |  N6  |			*/
			/*	|______|______|______	*/
			/*  |	   |	  |		 |	*/
			/*	|  I14 |  N7  |	 **  |	*/
			/*	|______|______|______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[3];				
			}
		}
	}

	//	Assigns the interactions for child1 of a box
	void assign_Child1_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+1;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N6  |  N5  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N7  |  **  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[7]	=	nC-1;
			tree[nL][nC].neighborNumbers[5]	=	nC+1;
			tree[nL][nC].neighborNumbers[6]	=	nC+2;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/************************************/
			/*				   		  ______	*/
			/*				  	     |		|	*/
			/*				         |	**  |	*/
			/*	 _____________       |______|  	*/
			/*	|	   |	  |					*/
			/*	|  O22 |  I15 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O23 |  I0  |					*/
			/*	|______|______|					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[23]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[22]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 		______		  	*/
			/*		   |	  |			*/
			/*	       |  **  |			*/
			/*	 ______|______|			*/
			/*	|	   |	  |			*/
			/*	|  N0  |  N1  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I1  |  I2  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's second neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |	  			*/
			/*	|______|_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  N2  |  I5  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*	       |  I3  |  I4  |	*/
			/*		   |______|______|	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's third neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  N4  |	 I7	 |  */
			/*	 ______|______|______|	*/
			/*	|	   |	  |	     |	*/
			/*	|  **  |  N3  |  I6  |  */
			/*	|______|______|______| 	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN != -1){
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  O14 |  O13 |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  I9  |  I8  |	*/
			/*		   |______|______|	*/
			/*				  			*/
			/*				  			*/
			/*	 ______					*/
			/*	|	   |				*/
			/*  |  **  |				*/
			/*  |______|				*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[13]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[14]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  O16 |  O15 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  I11 |  I10 |		*/
			/*	|______|______|		*/
			/*						*/
			/*						*/
			/*		    ______		*/
			/* 		   |	  |		*/
			/*		   |  **  |		*/
			/*		   |______|		*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[15]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[16]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/************************************/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O18 |  O17 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O19 |  I12 |					*/
			/*	|______|______|					*/
			/*									*/
			/*									*/
			/*				   		  ______	*/
			/*  			  		 |		|	*/
			/*				  		 |	** 	|	*/
			/*				  		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[19]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[17]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[18]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/************************************/
			/*									*/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O20 |  I13 |					*/
			/*	|______|______|		  ______	*/
			/*  |	   |	  |		 |		|	*/
			/*	|  O21 |  I14 |	 	 |	**  |	*/
			/*	|______|______|		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[21]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[20]	=	tree[j][nN].childrenNumbers[3];				
			}
		}
	}

	//	Assigns the interactions for child2 of a box
	void assign_Child2_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+2;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N7  |  **  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N0  |  N1  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[0]	=	nC-2;
			tree[nL][nC].neighborNumbers[1]	=	nC-1;
			tree[nL][nC].neighborNumbers[7]	=	nC+1;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/************************************/
			/*				   		  ______	*/
			/*				  	     |		|	*/
			/*				         |	**  |	*/
			/*				         |______|  	*/
			/*									*/
			/*									*/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O23 |  I0  |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O0  |  O1  |					*/
			/*	|______|______|					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[23]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 		______		  	*/
			/*		   |	  |			*/
			/*	       |  **  |			*/
			/*	 	   |______|			*/
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I1  |  I2  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O2  |  O3  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's second neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |	  			*/
			/*	|______|				*/
			/*							*/
			/*							*/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  I3  |  I4  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*	       |  O4  |  O5  |	*/
			/*		   |______|______|	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[5]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's third neighbor
		{
			/****************************/
			/*	 ____________________	*/
			/*	|	   |	  |	     |	*/
			/*	|  **  |  N3  |	 I6	 |  */
			/*	|______|______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  N2  |  I5  |  */
			/*		   |______|______| 	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  I9  |  I8  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  N4  |  I7  |	*/
			/*	 ______|______|______|	*/
			/*	|	   |	  			*/
			/*	|  **  |	  			*/
			/*	|______|	  			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  I11 |  I10 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N6  |  N5  |		*/
			/*	|______|______|		*/
			/*		   |	  |		*/
			/*		   |  **  |		*/
			/*		   |______|		*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/************************************/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O19 |  I12 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O20 |  I13 |					*/
			/*	|______|______|		  ______	*/
			/*  			  		 |		|	*/
			/*				  		 |	** 	|	*/
			/*				  		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[20]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[19]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/************************************/
			/*									*/
			/*	 _____________		  ______	*/
			/*	|	   |	  |		 |	    |	*/
			/*	|  O21 |  I14 |		 |	**	|	*/
			/*	|______|______|		 |______|	*/
			/*  |	   |	  |		 			*/
			/*	|  O22 |  I15 |	 	 			*/
			/*	|______|______|		 			*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[22]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[21]	=	tree[j][nN].childrenNumbers[3];				
			}
		}
	}

	//	Assigns the interactions for child3 of a box
	void assign_Child3_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+3;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  **  |  N3  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N1  |  N2  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[1]	=	nC-3;
			tree[nL][nC].neighborNumbers[2]	=	nC-2;
			tree[nL][nC].neighborNumbers[3]	=	nC-1;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/****************************/
			/*				   ______	*/
			/*				  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|  */
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I0  |  I1  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O1  |  O2  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |				*/
			/*	|______|				*/
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I2  |  I3  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O3  |  O4  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's second neighbor
		{
			/************************************/
			/*	 ______		  					*/
			/*	|	   |						*/
			/*	|  **  |	  					*/
			/*	|______|						*/
			/*									*/
			/*									*/
			/*				   _____________	*/
			/*		   		  |	     |	    |	*/
			/*		   		  |  I4  |  O7  |	*/
			/*		   		  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*	       		  |  O5  |  O6  |	*/
			/*		   		  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's third neighbor
		{
			/************************************/
			/*	 ______		   _____________	*/
			/*	|	   |	  |	     |		|	*/
			/*	|  **  |      |	 I6	 |  O9	|	*/
			/*	|______|	  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I5  |  O8  |  	*/
			/*		   		  |______|______| 	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[8]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[9]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I8  |  O11 |	*/
			/*		   		  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I7  |  O10 |	*/
			/*	 ______	      |______|______|	*/
			/*	|	   |	  					*/
			/*	|  **  |	  					*/
			/*	|______|	  					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[11]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  I10 |  I9  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N5  |  N4  |		*/
			/*	|______|______|		*/
			/*	|	   |			*/
			/*	|  **  |			*/
			/*	|______|			*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/****************************/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I12 |  I11 |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I13 |  N6  |			*/
			/*	|______|______|______	*/
			/*  			  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/****************************/
			/*							*/
			/*	 ____________________	*/
			/*	|	   |	  |		 |	*/
			/*	|  I14 |  N7  |	 **	 |	*/
			/*	|______|______|______|	*/
			/*  |	   |	  |		 	*/
			/*	|  I15 |  N0  |	 	 	*/
			/*	|______|______|		 	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[3];				
			}
		}
	}

	//	Assigns the interactions for the children of a box
	void assign_Box_Interactions(int j, int k) {
		assign_Child0_Interaction(j,k);
		assign_Child1_Interaction(j,k);
		assign_Child2_Interaction(j,k);
		assign_Child3_Interaction(j,k);
	}

	//	Assigns the interactions for the children all boxes at a given level
	void assign_Level_Interactions(int j) {
		#pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			assign_Box_Interactions(j,k);
		}
	}

	//	Assigns the interactions for the children all boxes in the tree
	void assign_Tree_Interactions() {
		for (int j=0; j<nLevels; ++j) {
			assign_Level_Interactions(j);
		}
	}

	//	Obtain the desired matrix
	void obtain_Desired_Operator(std::vector<pts2D>& shiftedChebNodes, Eigen::MatrixXd& T) {
		T	=	Eigen::MatrixXd(rank,rank);
		for (int i=0; i<rank; ++i) {
			for (int j=0; j<rank; ++j) {
				T(i,j)	=	K->getInteraction(standardChebNodes[i], shiftedChebNodes[j], a);
			}
		}
	}

	void obtain_Desired_Leaf_Operator(std::vector<pts2D>& shiftedChebNodes, Eigen::MatrixXd& T) {
		T	=	Eigen::MatrixXd(rank,rank);
		for (int i=0; i<rank; ++i) {
			for (int j=0; j<rank; ++j) {
				T(i,j)	=	K->getInteraction(leafChebNodes[i], shiftedChebNodes[j], a);
			}
		}
	}

	//	Assemble FMM Operators
	void assemble_Operators_FMM() {
		std::vector<pts2D> shiftedChebNodes;
		//	Assemble Outer Interactions
		for (int l=0; l<6; ++l) {
			shiftedChebNodes	=	shift_Cheb_Nodes(l-3,-3);
			obtain_Desired_Operator(shiftedChebNodes, M2LOuter[l]);
		}
		for (int l=0; l<6; ++l) {
			shiftedChebNodes	=	shift_Cheb_Nodes(3,l-3);
			obtain_Desired_Operator(shiftedChebNodes, M2LOuter[l+6]);
		}
		for (int l=0; l<6; ++l) {
			shiftedChebNodes	=	shift_Cheb_Nodes(3-l,3);
			obtain_Desired_Operator(shiftedChebNodes, M2LOuter[l+12]);
		}
		for (int l=0; l<6; ++l) {
			shiftedChebNodes	=	shift_Cheb_Nodes(-3,3-l);
			obtain_Desired_Operator(shiftedChebNodes, M2LOuter[l+18]);
		}
		//	Assemble Inner Interactions
		for (int l=0; l<4; ++l) {
			shiftedChebNodes	=	shift_Cheb_Nodes(l-2,-2);
			obtain_Desired_Operator(shiftedChebNodes, M2LInner[l]);
		}
		for (int l=0; l<4; ++l) {
			shiftedChebNodes	=	shift_Cheb_Nodes(2,l-2);
			obtain_Desired_Operator(shiftedChebNodes, M2LInner[l+4]);
		}
		for (int l=0; l<4; ++l) {
			shiftedChebNodes	=	shift_Cheb_Nodes(2-l,2);
			obtain_Desired_Operator(shiftedChebNodes, M2LInner[l+8]);
		}
		for (int l=0; l<4; ++l) {
			shiftedChebNodes	=	shift_Cheb_Nodes(-2,2-l);
			obtain_Desired_Operator(shiftedChebNodes, M2LInner[l+12]);
		}

		//	Assigning Leaf Interactions
		for (int k=0; k<rank; ++k) {
			pts2D temp;
			temp.x	=	boxRadius[nLevels]*standardChebNodes[k].x;
			temp.y	=	boxRadius[nLevels]*standardChebNodes[k].y;
			leafChebNodes.push_back(temp);
		}
		//	Assemble Neighbor Interactions
		double neighborDistance	=	2.0*boxRadius[nLevels];
		{
			shiftedChebNodes	=	shift_Leaf_Cheb_Nodes(-neighborDistance,-neighborDistance);
			obtain_Desired_Leaf_Operator(shiftedChebNodes, neighborInteraction[0]);
		}
		{
			shiftedChebNodes	=	shift_Leaf_Cheb_Nodes(0,-neighborDistance);
			obtain_Desired_Leaf_Operator(shiftedChebNodes, neighborInteraction[1]);
		}
		{
			shiftedChebNodes	=	shift_Leaf_Cheb_Nodes(neighborDistance,-neighborDistance);
			obtain_Desired_Leaf_Operator(shiftedChebNodes, neighborInteraction[2]);
		}
		{
			shiftedChebNodes	=	shift_Leaf_Cheb_Nodes(neighborDistance,0);
			obtain_Desired_Leaf_Operator(shiftedChebNodes, neighborInteraction[3]);
		}
		{
			shiftedChebNodes	=	shift_Leaf_Cheb_Nodes(neighborDistance,neighborDistance);
			obtain_Desired_Leaf_Operator(shiftedChebNodes, neighborInteraction[4]);
		}
		{
			shiftedChebNodes	=	shift_Leaf_Cheb_Nodes(0,neighborDistance);
			obtain_Desired_Leaf_Operator(shiftedChebNodes, neighborInteraction[5]);
		}
		{
			shiftedChebNodes	=	shift_Leaf_Cheb_Nodes(-neighborDistance,neighborDistance);
			obtain_Desired_Leaf_Operator(shiftedChebNodes, neighborInteraction[6]);
		}
		{
			shiftedChebNodes	=	shift_Leaf_Cheb_Nodes(-neighborDistance,0);
			obtain_Desired_Leaf_Operator(shiftedChebNodes, neighborInteraction[7]);
		}
		//	Assemble Self Interactions
		{
			obtain_Desired_Leaf_Operator(leafChebNodes, selfInteraction);
		}
	}
	
	void assign_Center_Location() {
		int J, K;
		tree[0][0].center.x	=	0.0;
		tree[0][0].center.y	=	0.0;
		for (int j=0; j<nLevels; ++j) {
			J	=	j+1;
			double shift	=	0.5*boxRadius[j];
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				K	=	4*k;
				tree[J][K].center.x		=	tree[j][k].center.x-shift;
				tree[J][K+1].center.x	=	tree[j][k].center.x+shift;
				tree[J][K+2].center.x	=	tree[j][k].center.x+shift;
				tree[J][K+3].center.x	=	tree[j][k].center.x-shift;

				tree[J][K].center.y		=	tree[j][k].center.y-shift;
				tree[J][K+1].center.y	=	tree[j][k].center.y-shift;
				tree[J][K+2].center.y	=	tree[j][k].center.y+shift;
				tree[J][K+3].center.y	=	tree[j][k].center.y+shift;
			}
		}
		#pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
			tree[nLevels][k].chebNodes	=	shift_Leaf_Cheb_Nodes(tree[nLevels][k].center.x, tree[nLevels][k].center.y);
		}
	}

	void assign_Leaf_Charges() {
		#pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
			tree[nLevels][k].multipoles	=	0.5*(Eigen::VectorXd::Ones(rank)+Eigen::VectorXd::Random(rank));
		}
	}

	void evaluate_All_M2M() {
		for (int j=nLevels-1; j>1; --j) {
			int J	=	j+1;
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				int K	=	4*k;
				tree[j][k].multipoles	=	M2M[0]*tree[J][K].multipoles+M2M[1]*tree[J][K+1].multipoles+M2M[2]*tree[J][K+2].multipoles+M2M[3]*tree[J][K+3].multipoles;
			}
		}
	}

	void evaluate_All_M2L() {
		for (int j=2; j<=nLevels; ++j) {
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				tree[j][k].locals	=	Eigen::VectorXd::Zero(rank);
				#ifdef HOMOG
					//	Inner well-separated clusters
					for (int l=0; l<16; ++l) {
						int nInner	=	tree[j][k].innerNumbers[l];
						if (nInner>-1) {
							tree[j][k].locals+=M2LInner[l]*tree[j][nInner].multipoles;
						}
					}
					//	Outer well-separated clusters
					for (int l=0; l<24; ++l) {
						int nOuter	=	tree[j][k].outerNumbers[l];
						if (nOuter>-1) {
							tree[j][k].locals+=M2LOuter[l]*tree[j][nOuter].multipoles;
						}
					}
					tree[j][k].locals*=boxHomogRadius[j];					
				#elif LOGHOMOG
					//	Inner well-separated clusters
					for (int l=0; l<16; ++l) {
						int nInner	=	tree[j][k].innerNumbers[l];
						if (nInner>-1) {
							tree[j][k].locals+=M2LInner[l]*tree[j][nInner].multipoles;
							tree[j][k].locals+=boxLogHomogRadius[j]*tree[j][nInner].multipoles.sum()*Eigen::VectorXd::Ones(rank);
						}
					}
					//	Outer well-separated clusters
					for (int l=0; l<24; ++l) {
						int nOuter	=	tree[j][k].outerNumbers[l];
						if (nOuter>-1) {
							tree[j][k].locals+=M2LOuter[l]*tree[j][nOuter].multipoles;
							tree[j][k].locals+=boxLogHomogRadius[j]*tree[j][nOuter].multipoles.sum()*Eigen::VectorXd::Ones(rank);
						}
					}
				#endif
			}
		}
	}

	void evaluate_All_L2L() {
		for (int j=2; j<nLevels; ++j) {
			int J	=	j+1;
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				int K	=	4*k;
				tree[J][K].locals+=L2L[0]*tree[j][k].locals;
				tree[J][K+1].locals+=L2L[1]*tree[j][k].locals;
				tree[J][K+2].locals+=L2L[2]*tree[j][k].locals;
				tree[J][K+3].locals+=L2L[3]*tree[j][k].locals;
			}
		}
	}

	void evaluate_Leaf() {
		if (nLevels <2) {
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
				tree[nLevels][k].locals	=	Eigen::VectorXd::Zero(rank);
			}
		}
		#pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
			for (int l=0; l<8; ++l) {
				int nNeighbor	=	tree[nLevels][k].neighborNumbers[l];
				if (nNeighbor > -1) {
					tree[nLevels][k].locals+=neighborInteraction[l]*tree[nLevels][nNeighbor].multipoles;
				}
			}
			tree[nLevels][k].locals+=selfInteraction*tree[nLevels][k].multipoles;
		}
	}

	double perform_Error_Check(int nBox) {
		Eigen::VectorXd potential	=	Eigen::VectorXd::Zero(rank);
		for (int l1=0; l1<rank; ++l1) {
			for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
				for (int l2=0; l2<rank; ++l2) {
					potential(l1)+=K->getInteraction(tree[nLevels][nBox].chebNodes[l1], tree[nLevels][k].chebNodes[l2], a)*tree[nLevels][k].multipoles(l2);
				}
			}
		}
		Eigen::VectorXd error(rank);
		for (int k=0; k<rank; ++k) {
			error(k)	=	fabs((potential-tree[nLevels][nBox].locals)(k)/potential(k));
		}
		return error.maxCoeff();
	}
};

#endif

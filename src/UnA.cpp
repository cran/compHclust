extern "C" {

#include <R.h>
	
	/* Private function prototypes. */
	void UnA(int* mMat, double* hts, int* nC, double* sMat, double* aMat);
	void RowDesc(int* mMat, int row, int nC, int* pnd, int* d);
    void EntDesc(int* mMat, int ent, int nC, int* pnd, int* d);

    /* UnA() returns an unnormalized version of the A matrix.
	Arguments:
	mMat - The merge matrix.  (A vector, columnwise)
	hts - The vector of delta h's.
	nC - The number of columns in the data.
	sMat - The step matrix of updates.  Initialized to be the identity.  (A vector, columnwise)
	aMat - The unnormalized A matrix.  Initialized to be (delta h1)*identity.
	Thus the R wrapper function needs to provide the arguments appropriately, and also divide
	the A matrix by the sum of the delta h's. */

	void UnA(int* mMat, double* hts, int* nC, double* sMat, double* aMat) {
		for (int i = 0; i < (*nC-2); i++) {
			int d[*nC];
			int nd = 0;
			int* pnd = &nd;
			RowDesc(mMat, i+1, *nC, pnd, d);
			for (int j = 0; j < nd; j++) {
				for (int k = 0; k < nd; k++) {
					sMat[d[j]+(*nC)*(d[k]-1)-1] = (1.0)/nd;
				}
			}
			for (int l = 0; l < ((*nC)*(*nC)); l++) {
				aMat[l] += hts[i+1]*(sMat[l]);
			}
		}
	}

	/* RowDesc() is a recursive function that determines the descendents of a row of the
	merge matrix.
	Arguments:
	mMat - The merge matrix.  (A vector, columnwise)
	row - A particular row of the merge matrix.
	nC - The number of columns in the data.
	pnd - A pointer to the number of descendents.
	d - The vector of descendents. */

	void RowDesc(int* mMat, int row, int nC, int* pnd, int* d) {
		EntDesc(mMat, row, nC, pnd, d);
		EntDesc(mMat, row+nC-1, nC, pnd, d);
	}

	/* EntDesc() is a recursive function that determines the descendents of an entry of the
	merge matrix.
	Arguments:
	mMat - The merge matrix.  (A vector, columnwise)
	ent - A particular entry of the merge matrix.
	nC - The number of columns in the data.
	pnd - A pointer to the number of descendents.
	d - The vector of descendents. */

    void EntDesc(int* mMat, int ent, int nC, int* pnd, int* d) {
        if (mMat[ent-1]<0) {
            d[*pnd] = -mMat[ent-1];
            (*pnd) += 1;
		} else {
            RowDesc(mMat, mMat[ent-1], nC, pnd, d);
		}
	}

}

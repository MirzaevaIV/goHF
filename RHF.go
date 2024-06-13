// RHF.go --  This file is part of goHF project.
// Mirzaeva Irina, 2023
//
//	goHF is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty
//	of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//	See the GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see http://www.gnu.org/licenses/
//
// ------------------------------------------------
package main

import "C"

import (
	_ "flag"
	"fmt"
	_ "log"
	"math"
	"time"
	_ "os"
	"runtime"
	"sync"
	_ "reflect"
	_ "strings"
	_ "strconv"

	_ "golang.org/x/exp/slices"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
)

type RHF struct {
	Occupied                int
	S, T, Ven               [][]float64
	Vnn                     float64
	VeeIdx                  []int
	VeeVal                  []float64
	S2Inv, H1, Cij, DensMat, G [][]float64
	F_list, DIIS_R []*mat.Dense
}

type MatrixData struct {
	i, j int
	data []float64
}

func (m *Molecule) RHFinit() RHF {
	tstart := time.Now()
	var result RHF

	result.Occupied = m.getNelec() / 2

	fmt.Println("Initializing RHF structure...")
	result.S, result.T, result.Ven = m.CalculateIntegrals_1e()
	tstop := time.Now()
	fmt.Println("------", "Time for building CINT data and 1e integrals calculation:", tstop.Sub(tstart))
	InfoLogger.Println("1e integrals done...",  tstop.Sub(tstart))
	tstart = tstop

	result.VeeIdx, result.VeeVal = m.ElecElec1()
	tstop = time.Now()
	fmt.Println("------", "Time for 2e integrals calculation:", tstop.Sub(tstart))
	InfoLogger.Println("2e integrals done...",  tstop.Sub(tstart))
	tstart = tstop

	result.S2Inv = MatrixSqrtInverse(result.S)

	result.BuildInitialGuess()
	result.BuildDensMat()

	result.Vnn = m.NucNuc()
	fmt.Println("Finished building DensMat:")
	return result
}

func (rhf *RHF) BuildInitialGuess() {
	n_basis := len(rhf.T)
	H1 := mat.NewDense(n_basis, n_basis, flatten(rhf.T))
	H1.Add(H1, mat.NewDense(n_basis, n_basis, flatten(rhf.Ven)))

	rhf.H1 = make([][]float64, n_basis)
	for i := range rhf.H1 {
		rhf.H1[i] = make([]float64, n_basis)
		copy(rhf.H1[i], H1.RawRowView(i))
	}

	SSqrtInv := mat.NewDense(n_basis, n_basis, flatten(rhf.S2Inv))

	H1.Mul(SSqrtInv, H1)
	H1.Mul(H1, SSqrtInv)

	H1Sym := mat.NewSymDense(n_basis, H1.RawMatrix().Data)
	var eigsym mat.EigenSym
	ok := eigsym.Factorize(H1Sym, true)
	if !ok {
		fmt.Println("Transformed H1 eigendecomposition failed")
	}

	var ev mat.Dense
	eigsym.VectorsTo(&ev)

	ev.Mul(SSqrtInv, &ev)

	vals := make([]float64, n_basis)
	eigsym.Values(vals)
	fmt.Println(vals)

	rhf.Cij = make([][]float64, n_basis)
	for i := range rhf.Cij {
		rhf.Cij[i] = ev.RawRowView(i)
	}

}

func (rhf *RHF) BuildDensMat() {
	if len(rhf.Cij) == 0{
		fmt.Println("Cannot build Density matrix. No Cij.")
		//break
	}
	n_basis := len(rhf.Cij)
	rhf.DensMat = make([][]float64, n_basis)
	for i := range rhf.DensMat {
		rhf.DensMat[i] = make([]float64, n_basis)
	}
	occ := 1.0
	for i := 0; i < n_basis; i++ {
		for j := 0; j < n_basis; j++ {
			for oo := 0; oo < rhf.Occupied; oo++ {
				rhf.DensMat[i][j] += occ * rhf.Cij[i][oo] * rhf.Cij[j][oo]
			}
		}
	}

}

func (rhf *RHF) SetGZeros(){
	n_basis := len(rhf.T)

	rhf.G = make([][]float64, n_basis)
	for i := range rhf.G {
		rhf.G[i] = make([]float64, n_basis)
	}
}

func (rhf *RHF) ProcessVeeVal(idx int, G [][]float64) {
	i, j, k, l := rhf.GetVeeIdxs(rhf.VeeIdx[idx])
	G[i][j] += 2*rhf.DensMat[k][l] * rhf.VeeVal[idx]
	G[i][k] -= rhf.DensMat[j][l] * rhf.VeeVal[idx]
	if i != j {
		G[j][i] += 2*rhf.DensMat[k][l] * rhf.VeeVal[idx]
		G[j][k] -= rhf.DensMat[i][l] * rhf.VeeVal[idx]
		if k!=l {
			G[j][i] += 2*rhf.DensMat[l][k] * rhf.VeeVal[idx]
			G[j][l] -= rhf.DensMat[i][k] * rhf.VeeVal[idx]
		}
	}
	if k!=l{
		G[i][j] += 2*rhf.DensMat[l][k] * rhf.VeeVal[idx]
		G[i][l] -= rhf.DensMat[j][k] * rhf.VeeVal[idx]
	}
}

func (rhf *RHF) ProcessVeeVal1(idx int, result *mat.Dense){
	i, j, k, l := rhf.GetVeeIdxs(rhf.VeeIdx[idx])
	result.Set(i, j, result.At(i, j)+2*rhf.DensMat[k][l] * rhf.VeeVal[idx])
	result.Set(i, k, result.At(i, k)-rhf.DensMat[j][l] * rhf.VeeVal[idx])
	if i != j {
		result.Set(j, i, result.At(j, i)+2*rhf.DensMat[k][l] * rhf.VeeVal[idx])
	    result.Set(j, k, result.At(j, k)-rhf.DensMat[i][l] * rhf.VeeVal[idx])
		if k!=l {
			result.Set(j, i, result.At(j, i)+2*rhf.DensMat[l][k] * rhf.VeeVal[idx])
	        result.Set(j, l, result.At(j, l)-rhf.DensMat[i][k] * rhf.VeeVal[idx])
		}
	}
	if k!=l{
		result.Set(i, j, result.At(i, j)+2*rhf.DensMat[l][k] * rhf.VeeVal[idx])
	    result.Set(i, l, result.At(i, l)-rhf.DensMat[j][k] * rhf.VeeVal[idx])
	}
}

func (rhf *RHF) ProcessVeeIdx(idxList []int, start int, G [][]float64) {
	for idx := range idxList {
		rhf.ProcessVeeVal((idx+start), G)
	}
}

func (rhf *RHF) ProcessVeeIdx1(idxList []int, start int, G *mat.Dense) {
	for idx := range idxList {
		rhf.ProcessVeeVal1((idx+start), G)
	}
}

func (rhf *RHF) BuildG() {
	if len(rhf.DensMat) == 0{
		rhf.BuildDensMat()
	}
	rhf.SetGZeros()
	n_basis := len(rhf.T)
	NVee := len(rhf.VeeIdx)
	maxGoroutines := runtime.GOMAXPROCS(-1)
	fmt.Println("NVee, maxGoroutines", NVee, maxGoroutines)
	var wg sync.WaitGroup
	if maxGoroutines>1{
		ListSize := NVee / (maxGoroutines)
		fmt.Println("ListSize", ListSize)
		Gparts := make([]*mat.Dense, maxGoroutines)
		for i := range Gparts {
			Gparts[i] = mat.NewDense(n_basis, n_basis, nil)
		}
		for j:=0; j<(maxGoroutines-1);j++{
			wg.Add(1)
		    go func(){
				defer wg.Done()
                rhf.ProcessVeeIdx1(rhf.VeeIdx[(j*ListSize):((j+1)*ListSize)], j*ListSize, Gparts[j])
			}()
	    }
		wg.Add(1)
		go func(){
			defer wg.Done()
		    rhf.ProcessVeeIdx1(rhf.VeeIdx[((maxGoroutines-1)*ListSize):], (maxGoroutines-1)*ListSize, Gparts[maxGoroutines-1])
		}()
		wg.Wait()
		result:= mat.NewDense(n_basis, n_basis, nil)
		for i := range Gparts {
			result.Add(result, Gparts[i])
		}
		for i := range rhf.G {
			rhf.G[i] = result.RawRowView(i)
		}
	} else{
		rhf.ProcessVeeIdx(rhf.VeeIdx, 0, rhf.G)
	} 
}


func (rhf *RHF) CalcEnergy() float64 {
	n_basis := len(rhf.T)
	res := 0.0
	for i := 0; i < n_basis; i++ {
		for j := 0; j < n_basis; j++ {
			res += rhf.DensMat[i][j] * (2*rhf.H1[i][j] + rhf.G[i][j])
		}
	}
	return res + rhf.Vnn
}

func (rhf *RHF) GetVeeIdxs(IdxVal int) (int, int, int, int){
	n_basis := len(rhf.T)
	i:= IdxVal/(n_basis*n_basis*n_basis)
	IdxVal = IdxVal%(n_basis*n_basis*n_basis)
	j:= IdxVal/(n_basis*n_basis)
	IdxVal = IdxVal%(n_basis*n_basis)
	k:= IdxVal/n_basis
	l := IdxVal%n_basis
	return i, j, k, l
}


func (rhf *RHF) BuildDIIS_R(F, S2inv *mat.Dense) {
	n_basis := len(rhf.T)
	term1 := mat.NewDense(n_basis, n_basis, nil)
	term2 := mat.NewDense(n_basis, n_basis, nil)
	S := mat.NewDense(n_basis, n_basis, flatten(rhf.S))
	DM := mat.NewDense(n_basis, n_basis, flatten(rhf.DensMat))
	term1.Mul(F, DM)
	term1.Mul(term1, S)
	term2.Mul(S, DM)
	term2.Mul(term2, F)
	term1.Sub(term1, term2)
	term1.Mul(S2inv, term1)
	term1.Mul(term1, S2inv)
	rhf.DIIS_R = append(rhf.DIIS_R, term1) //DIIS Residual: diis_r = A.dot(F.dot(D).dot(S) - S.dot(D).dot(F)).dot(A)
}

func (rhf *RHF) CalcdRMS() float64 {
	result := 0.0
	res := mat.DenseCopyOf(rhf.DIIS_R[len(rhf.DIIS_R)-1])
	res.MulElem(res, res)
	result = math.Sqrt(stat.Mean(res.RawMatrix().Data, nil))
	return result
}

func (rhf *RHF) BuildB() *mat.Dense { //https://github.com/psi4/psi4numpy/blob/master/Tutorials/03_Hartree-Fock/3b_rhf-diis.ipynb
	B_dim := len(rhf.F_list) + 1
	r_dim := len(rhf.T)
	result := mat.NewDense(B_dim, B_dim, nil)

	for i:=0; i<(B_dim-1); i++{
		result.Set(i, B_dim-1, -1)
		result.Set(B_dim-1, i, -1)
	}

	for i := range rhf.F_list{
		for j := range rhf.F_list{
			b := mat.NewDense(r_dim, r_dim, nil)
			b.MulElem(rhf.DIIS_R[i], rhf.DIIS_R[j])
			result.Set(i, j, mat.Sum(b))
		}
	}
	return result
}


func (rhf *RHF) SCF_DIIS() float64 {
	tstart := time.Now()
	TolE := 1e-6    //Tolerance
	TolD := 1e-3
	MaxSteps := 20 //50 //SCF Maximum Steps
	n_basis := len(rhf.H1)
	res := 0.0
	E_prev := 0.0
	dRMS := 0.0

	H1 := mat.NewDense(n_basis, n_basis, flatten(rhf.H1))
	SSqrtInv := mat.NewDense(n_basis, n_basis, flatten(rhf.S2Inv))

	tstop := time.Now()
	//fmt.Println("***---***---", "Time for initialization of SCF procedure:", tstop.Sub(tstart))
	tstart = tstop

	for i := 0; i < MaxSteps; i++ {
		E_prev = res
		rhf.BuildG()

		res = rhf.CalcEnergy()

		F := mat.NewDense(n_basis, n_basis, flatten(rhf.G))
		F.Add(F, H1)
			
		rhf.F_list = append(rhf.F_list, mat.DenseCopyOf(F))
		rhf.BuildDIIS_R(F, SSqrtInv)
    	dRMS = rhf.CalcdRMS()

		OutputLogger.Println("Iteration ", i+1, ". Energy = ", res, ", dE = ", E_prev-res, ", dRMS = ", dRMS)
		fmt.Println("Iteration ", i+1, ". Energy = ", res, ", dE = ", E_prev-res, ", dRMS = ", dRMS)
		if (math.Abs(E_prev-res) < TolE) && (dRMS < TolD) {
			OutputLogger.Println("SCF converged after step ", i+1)
			fmt.Println("SCF converged after step ", i+1)
			return res
		}

		if i>0{
			bmat := rhf.BuildB()

			rhs := mat.NewVecDense((len(rhf.F_list) + 1), nil)
			rhs.SetVec(len(rhf.F_list) , -1)

			var lu mat.LU
		    lu.Factorize(bmat)
		    var coefs mat.VecDense
		    if err := lu.SolveVecTo(&coefs, false, rhs); err != nil {
			    continue
		    }
			F = mat.NewDense(n_basis, n_basis, nil)
			for j := range rhf.F_list{
				fpart := mat.NewDense(n_basis, n_basis, nil)
				fpart.Scale(coefs.AtVec(j), rhf.F_list[j])
				F.Add(F, fpart)
			}
		} 

		F.Mul(F, SSqrtInv)
		F.Mul(SSqrtInv, F)
		FSym := mat.NewSymDense(n_basis, F.RawMatrix().Data)
		var eigsym mat.EigenSym
		var ev mat.Dense
		ok := eigsym.Factorize(FSym, true)
		if !ok {
			fmt.Println("Transformed F eigendecomposition failed")
		}

		eigsym.VectorsTo(&ev)
		ev.Mul(SSqrtInv, &ev)

		for i := range rhf.Cij {
			rhf.Cij[i] = ev.RawRowView(i)
		}

		rhf.BuildDensMat()

		tstop = time.Now()
		fmt.Println("***---***---", "Time for SCF with DIIS iteration #", i+1, ":",  tstop.Sub(tstart))
		tstart = tstop
	}

	OutputLogger.Println("Warning! SCF NOT converged after step ", MaxSteps)
	fmt.Println("Warning! SCF NOT converged after step ", MaxSteps)
	return res
}


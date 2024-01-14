// main.go --  This file is part of goHF project.
// Mirzaeva Irina, 2023
//
//  goHF is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty
//  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see http://www.gnu.org/licenses/
//------------------------------------------------
package main

import "C"

import (
	_ "flag"
	"fmt"
	"log"
	"math"
	"os"
	"strings"

	_ "golang.org/x/exp/slices"
	"gonum.org/v1/gonum/mat"
)

var (
	WarningLogger *log.Logger
	InfoLogger    *log.Logger
	ErrorLogger   *log.Logger
	OutputLogger  *log.Logger
)

var ElemData Mendeleev

var a_B = 0.52917720859

func init() {
	ElemData.build()
}

func initLog(fname string) {
	file, err := os.OpenFile(fname, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatal(err)
	}

	InfoLogger = log.New(file, "INFO: ", log.Ldate|log.Ltime)
	WarningLogger = log.New(file, "WARNING: ", log.Ldate|log.Ltime)
	ErrorLogger = log.New(file, "ERROR: ", log.Ldate|log.Ltime|log.Lshortfile)
	OutputLogger = log.New(file, "", 0)
}

func appInfo() {
	OutputLogger.Println("\n              __  __  ____      |\n             /\\ \\/\\ \\/\\  __\\    |" +
		" Author: Mirzaeva Irina Valerievna\n   __     ___\\ \\ \\_\\ \\ \\ \\_/    | email: dairdre@gmail.com\n" +
		" /'_ `\\  / __`\\ \\  _  \\ \\  _\\   | Nikolaev Institute of Inorganic Chemistry SB RAS" +
		" (http://niic.nsc.ru/)\n/\\ \\L\\ \\/\\ \\L\\ \\ \\ \\ \\ \\ \\ \\/   | Novosibirsk, Russia" +
		"\n\\ \\____ \\ \\____/\\ \\_\\ \\_\\ \\_\\   | HF stands for Himicheskaya Fizika\n \\/___L\\" +
		" \\/___/  \\/_/\\/_/\\/_/   | Have Fun!!!\n   /\\____/                      |\n   \\_/__/                       |\n\n")
}

func printOutputDelimiter() {
	OutputLogger.Println(strings.Repeat("-", 70))
}

func processInput(data []string) Molecule {
	var atoms, basis bool
	var atom_start, atom_end int
	var basisName string
	var mol Molecule
	for i := 0; i < len(data); i++ {
		words := strings.Fields(data[i])
		if len(words) > 0 {
			if strings.ToLower(words[0]) == "atoms" {
				atoms = true
				atom_start = i
				atom_end = findBlockEnd(i, data, "Atoms")
				OutputLogger.Print("Atoms block found at lines ", atom_start, " -- ", atom_end, ".")
			}
			if strings.ToLower(words[0]) == "basis" {
				basis = true
				basisName = data[i+1]
				_ = findBlockEnd(i, data, "Basis")
				OutputLogger.Print("Basis block found.", basisName)
			}
		}
	}
	if !atoms {
		ErrorLogger.Fatal("No Atoms found.")
	} else {
		mol.addAtoms(data, atom_start+1, atom_end-1)
	}
	if !basis {
		OutputLogger.Println("No Basis found. Using default basis: STO-3G.")
	} else {
		mol.getBasis(basisName)
	}
	return mol
}

func findBlockEnd(n int, data []string, bname string) int {
	for i := n; i < len(data); i++ {
		words := strings.Fields(data[i])
		if len(words) > 0 {
			if strings.ToLower(words[0]) == "end" {
				return i
			}
		}
	}
	ErrorLogger.Fatal("No end of block " + bname + ".")
	return 0
}

func BuildInitialGuess(S2Inv, T, Ven [][]float64) ([][]float64, [][]float64) {
	n_basis := len(T)
	H1 := mat.NewDense(n_basis, n_basis, flatten(T))
	H1.Add(H1, mat.NewDense(n_basis, n_basis, flatten(Ven)))

	resH1 := mat.NewDense(n_basis, n_basis, flatten(T))
	resH1.Copy(H1)

	SSqrtInv := mat.NewDense(n_basis, n_basis, flatten(S2Inv))

	H1.Mul(SSqrtInv, H1)
	H1.Mul(H1, SSqrtInv)

	H1Sym := mat.NewSymDense(n_basis, H1.RawMatrix().Data)
	var eigsym mat.EigenSym
	ok := eigsym.Factorize(H1Sym, true)
	if !ok {
		fmt.Println("Transformed H1 eigendecomposition failed")
	}

	var eig mat.Eigen
	ok = eig.Factorize(H1, mat.EigenRight)
	if !ok {
		fmt.Println("Transformed H1 eigendecomposition failed")
	}

	var ev1 mat.CDense
	eig.VectorsTo(&ev1)
	vals1 := make([]complex128, n_basis)
	eig.Values(vals1)
	vals1r := make([]float64, n_basis)
	for i := range vals1 {
		vals1r[i] = real(vals1[i])
	}
	fmt.Println(vals1r)
	vecr := make([][]float64, n_basis)
	for i := range vecr {
		vecr[i] = make([]float64, n_basis)
	}
	for i := range vecr {
		for j := range vecr[i] {
			vecr[i][j] = real(ev1.At(i, j))
		}
	}

	var ev mat.Dense
	eigsym.VectorsTo(&ev)

	ev.Mul(SSqrtInv, &ev)

	vals := make([]float64, n_basis)
	eigsym.Values(vals)

	C0 := make([][]float64, n_basis)
	for i := range C0 {
		C0[i] = ev.RawRowView(i)
	}

	rH1 := make([][]float64, n_basis)
	for i := range rH1 {
		rH1[i] = resH1.RawRowView(i)
	}
	return rH1, C0
}

func BuildDensMat(MOs [][]float64, nelec int) [][]float64 {
	n_basis := len(MOs)
	result := make([][]float64, n_basis)
	for i := range result {
		result[i] = make([]float64, n_basis)
	}
	occ := 1.0
	for i := 0; i < n_basis; i++ {
		for j := 0; j < n_basis; j++ {
			for oo := 0; oo < nelec; oo++ {
				C := MOs[i][oo]
				C_dag := MOs[j][oo]
				result[i][j] += occ * C * C_dag
			}
		}
	}
	return result
}

func CalcEnergy(H1, D, G [][]float64, Enuc float64) float64 {
	n_basis := len(H1)
	res := 0.0
	for i := 0; i < n_basis; i++ {
		for j := 0; j < n_basis; j++ {
			res += D[i][j] * (2*H1[i][j] + G[i][j])
		}
	}
	return res + Enuc
}

func BuildG(den_mat [][]float64, Vee [][][][]float64) [][]float64 {
	res := make([][]float64, len(den_mat))
	for i := range res {
		res[i] = make([]float64, len(den_mat))
	}
	for i := range Vee {
		for j := range Vee {
			for k := range Vee {
				for l := range Vee {
					J := Vee[i][j][k][l]
					K := Vee[i][k][j][l]
					res[i][j] += den_mat[k][l] * (2*J - K)
				}
			}
		}
	}
	return res
}

func BuildJ(den_mat [][]float64, Vee [][][][]float64) [][]float64 {
	res := make([][]float64, len(den_mat))
	for i := range res {
		res[i] = make([]float64, len(den_mat))
	}
	for i := range Vee {
		for j := range Vee {
			for k := range Vee {
				for l := range Vee {
					J := Vee[i][j][k][l]
					res[i][j] += den_mat[k][l] * (J)
				}
			}
		}
	}
	return res
}

func BuildK(den_mat [][]float64, Vee [][][][]float64) [][]float64 {
	res := make([][]float64, len(den_mat))
	for i := range res {
		res[i] = make([]float64, len(den_mat))
	}
	for i := range Vee {
		for j := range Vee {
			for k := range Vee {
				for l := range Vee {
					K := Vee[i][k][j][l]
					res[i][j] += den_mat[k][l] * (K)
				}
			}
		}
	}
	return res
}

func SCF(H1m, DM, S2inv [][]float64, Vee [][][][]float64, nelec int, Enuc float64) float64 {
	Tol := 1e-8    //Tolerance
	MaxSteps := 50 //SCF Maximum Steps
	n_basis := len(H1m)
	res := 0.0
	E_prev := 0.0

	H1 := mat.NewDense(n_basis, n_basis, flatten(H1m))
	SSqrtInv := mat.NewDense(n_basis, n_basis, flatten(S2inv))

	for i := 0; i < MaxSteps; i++ {
		E_prev = res
		G := BuildG(DM, Vee)

		res = CalcEnergy(H1m, DM, G, Enuc)
		OutputLogger.Println("Iteration ", i, ". Energy = ", res, ", dE = ", E_prev-res)
		fmt.Println("Iteration ", i, ". Energy = ", res, ", dE = ", E_prev-res)
		if math.Abs(E_prev-res) < Tol {
			OutputLogger.Println("SCF converged after step ", i)
			fmt.Println("SCF converged after step ", i)
			return res
		}

		F := mat.NewDense(n_basis, n_basis, flatten(G))

		F.Add(F, H1)

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

		MOs := make([][]float64, n_basis)
		for i := range MOs {
			MOs[i] = ev.RawRowView(i)
		}

		DM = BuildDensMat(MOs, nelec)
	}

	OutputLogger.Println("Warning! SCF NOT converged after step ", MaxSteps)
	fmt.Println("Warning! SCF NOT converged after step ", MaxSteps)
	return res
}

func main() {

	var inpFname, outFname string
	if len(os.Args) > 1 {
		inpFname = os.Args[1]
		split_inpFname := strings.Split(inpFname, ".")
		fExt := split_inpFname[len(split_inpFname)-1]
		outFname = inpFname[0:(len(inpFname)-len(fExt))] + "out"
		fmt.Println("Output file: ", outFname)
	} else {
		log.Fatal("No input file.")
	}

	initLog(outFname)

	InfoLogger.Println("Starting goHF...")
	appInfo()

	WarningLogger.Println("This is an experimental program on an early stage of development.")

	OutputLogger.Println("\n")
	OutputLogger.Println("Input file content:")
	printOutputDelimiter()
	inpData, err := ReadFileLines(inpFname)
	if err != nil {
		ErrorLogger.Println("Cannot read input file: ", err)
	}
	for _, i := range inpData {
		OutputLogger.Println(i)
	}
	printOutputDelimiter()

	mol := processInput(inpData)
	S, T, Ven, Vee := mol.CalculateIntegrals()

	S2Inv := MatrixSqrtInverse(S)

	H1, C0 := BuildInitialGuess(S2Inv, T, Ven)

	DensMat := BuildDensMat(C0, mol.getNelec()/2)

	G := BuildG(DensMat, Vee)
	E := CalcEnergy(H1, DensMat, G, mol.NucNuc())
	OutputLogger.Println("Initial energy = ", E, " a.u.")
	fmt.Println("Initial energy = ", E, " a.u.")

	EEE := SCF(H1, DensMat, S2Inv, Vee, mol.getNelec()/2, mol.NucNuc())

	OutputLogger.Println("Nuclei Repulsion Energy: ", mol.NucNuc(), " a.u.")
	printOutputDelimiter()
	fmt.Println("Nuc energy = ", mol.NucNuc(), " a.u.")

	fmt.Println("Final total energy = ", EEE, " a.u.")
	OutputLogger.Println("Final total energy = ", EEE, " a.u.")
	printOutputDelimiter()

	OutputLogger.Println("\n")
	InfoLogger.Println("Exiting goHF...")
	fmt.Println("goHF done.")
}

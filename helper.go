// helper.go --  This file is part of goHF project.
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

import (
	"C"
	"bufio"
	"fmt"
	"os"
	"runtime"
	_ "runtime/debug"
	"strconv"
)
import (
	"math"

	"gonum.org/v1/gonum/mat"
)

func ReadFileLines(fname string) ([]string, error) {
	var result []string
	var err error

	file, err := os.Open(fname)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		result = append(result, scanner.Text())
	}
	err = scanner.Err()

	return result, err
}

func TxtFileFrom2DSlice(data [][]float64, fname string) {
	var ftext string
	for i := 0; i < len(data); i++ {
		for j := 0; j < len(data[i]); j++ {
			ftext += fmt.Sprintf("%12.6f", data[i][j])
		}
		ftext += "\n"
	}
	err := os.WriteFile(fname, []byte(ftext), 0644)
	if err != nil {
		fmt.Println(err)
	}

}

func TxtFilesFrom4DSlice(data [][][][]float64, basename string) {
	os.Mkdir(basename, 644)
	for i := 0; i < len(data); i++ {
		os.Mkdir("./"+basename+"/"+basename+strconv.Itoa(i), 644)
		for j := 0; j < len(data[i]); j++ {
			fname := "./" + basename + "/" + basename + strconv.Itoa(i) + "/" + strconv.Itoa(j) + ".txt"
			TxtFileFrom2DSlice(data[i][j], fname)
		}
	}
}

func CConvertInt(a [][]int) [][]C.int {
	res := make([][]C.int, len(a))
	for i := range res {
		res[i] = make([]C.int, len(a[i]))
	}
	for i := range a {
		for j := range a[i] {
			res[i][j] = C.int(a[i][j])
		}
	}
	return res
}

func CConvertDouble(a []float64) []C.double {
	res := make([]C.double, len(a))
	for i := range res {
		res[i] = C.double(a[i])
	}
	return res
}

func flatten(arr [][]float64) []float64 {
	dim := len(arr)
	res := make([]float64, dim*dim)
	for i := range arr {
		for j := range arr[i] {
			res[i*dim+j] = arr[i][j]
		}
	}
	return res
}

func PrintMat(M [][]float64) {
	aaa := mat.NewDense(len(M), len(M), flatten(M))
	PrintDense(aaa)
}

func PrintDense(D *mat.Dense){
	fa := mat.Formatted(D, mat.Prefix("    "), mat.Squeeze())
	fmt.Printf("    %.8f\n", fa)
}

func MatrixSqrtInverse(S [][]float64) [][]float64 {
	n_basis := len(S)
	Smat := mat.NewSymDense(n_basis, flatten(S))
	var eigsym mat.EigenSym
	ok := eigsym.Factorize(Smat, true)
	if !ok {
		fmt.Println("S eigendecomposition failed")
	}
	var ev mat.Dense
	eigsym.VectorsTo(&ev)
	sqrtVec := make([]float64, n_basis)
	for i := range eigsym.Values(nil) {
		sqrtVec[i] = math.Sqrt(eigsym.Values(nil)[i])
	}
	diagM := mat.NewDiagDense(n_basis, sqrtVec)
	var SSqrtInv mat.Dense
	SSqrtInv.Mul(&ev, diagM)
	ev.Inverse(&ev)
	SSqrtInv.Mul(&SSqrtInv, &ev)
	SSqrtInv.Inverse(&SSqrtInv)

	result := make([][]float64, n_basis)
	for i := range result {
		result[i] = SSqrtInv.RawRowView(i)
	}
	return result
}

func MyMemDebug() {
	fmt.Println("-----------!!!!!!!!Enter MyMemDebug!!!!!!!!--------------")
	var memStats runtime.MemStats

	runtime.ReadMemStats(&memStats)

	fmt.Printf("Alloc: %d bytes\n", memStats.Alloc)
	fmt.Printf("TotalAlloc: %d bytes\n", memStats.TotalAlloc)
	fmt.Printf("HeapAlloc: %d bytes\n", memStats.HeapAlloc)
	fmt.Printf("HeapSys: %d bytes\n", memStats.HeapSys)
	fmt.Println("------------------------------------------!--------------")
}

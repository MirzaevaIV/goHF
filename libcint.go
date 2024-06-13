// libcint.go --  This file is part of goHF project.
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

/*
#cgo CFLAGS: -g -Wall
#cgo LDFLAGS: ./cfiles/libcint.a
#include <stdlib.h>
#include <math.h>
#include "./cfiles/cint.h"
#cgo LDFLAGS: -ltcmalloc -lm  -lquadmath
int cint1e_kin_sph(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env);
int cint1e_nuc_sph(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env);
int cint1e_ovlp_sph(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env);
int cint2e_sph(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt);
*/
import "C"
import (
	"fmt"
	"math"
	"runtime"
	"sync"
)

type CintData struct {
	di, dj, ptr_i, ptr_j int
	data                 []C.double
}

type Cint2Data struct {
	di, dj, dk, dl, ptr_i, ptr_j, ptr_k, ptr_l int
	data                                       []C.double
}

func (m *Molecule) CintOvlp(i, j, ptr_i, ptr_j int) CintData {
	di := m.Shells[i]
	dj := m.Shells[j]
	buff_size := di * dj
	shls := [2]C.int{C.int(i), C.int(j)}
	buf := make([]C.double, buff_size)

	_ = C.cint1e_ovlp_sph(&buf[0], &shls[0], &m.catm[0][0], m.natm, &m.cbas[0][0], m.nbas, &m.cenv[0])

	return CintData{di, dj, ptr_i, ptr_j, buf[:buff_size]}
}

func (m *Molecule) CintKinetic(i, j, ptr_i, ptr_j int) CintData {
	di := m.Shells[i]
	dj := m.Shells[j]
	buff_size := di * dj
	shls := [2]C.int{C.int(i), C.int(j)}
	buf := make([]C.double, buff_size)

	_ = C.cint1e_kin_sph(&buf[0], &shls[0], &m.catm[0][0], m.natm, &m.cbas[0][0], m.nbas, &m.cenv[0])
	return CintData{di, dj, ptr_i, ptr_j, buf[:buff_size]}
}

func (m *Molecule) CintElecNuc(i, j, ptr_i, ptr_j int) CintData {
	di := m.Shells[i]
	dj := m.Shells[j]
	buff_size := di * dj
	shls := [2]C.int{C.int(i), C.int(j)}
	buf := make([]C.double, buff_size)

	_ = C.cint1e_nuc_sph(&buf[0], &shls[0], &m.catm[0][0], m.natm, &m.cbas[0][0], m.nbas, &m.cenv[0])
	return CintData{di, dj, ptr_i, ptr_j, buf[:buff_size]}
}

func (m *Molecule) Ovlp() [][]float64 {
	var result [][]float64

	result = make([][]float64, m.nao)
	for i := range result {
		result[i] = make([]float64, m.nao)
	}

	chan1 := make(chan CintData)
	maxGoroutines := runtime.GOMAXPROCS(-1)
	guard := make(chan struct{}, maxGoroutines)

	go func() {
		res_ptri := 0
		for i, _ := range m.Shells {
			res_ptrj := 0
			for j := 0; j < i+1; j++ {
				go func(ii, jj, prti, ptrj int) {
					chan1 <- m.CintOvlp(ii, jj, prti, ptrj)
					<-guard
				}(i, j, res_ptri, res_ptrj)
				res_ptrj += m.Shells[j]
			}
			res_ptri += m.Shells[i]
		}
	}()

	var wg sync.WaitGroup
	wg.Add(1)

	go func() {
		defer wg.Done()
		res_ptri := 0
		for i, _ := range m.Shells {
			res_ptrj := 0
			for j := 0; j < i+1; j++ {
				data := <-chan1

				c := 0
				for dj := 0; dj < data.dj; dj++ {
					for di := 0; di < data.di; di++ {
						result[data.ptr_i+di][data.ptr_j+dj] = float64(data.data[c])
						result[data.ptr_j+dj][data.ptr_i+di] = float64(data.data[c])
						c++
					}
				}
				res_ptrj += m.Shells[j]
			}
			res_ptri += m.Shells[i]
		}
	}()
	wg.Wait()

	fmt.Println("Overlap integrals calculation done. Var 1.")
	return result
}

func (m *Molecule) Kinetic() [][]float64 {
	var result [][]float64

	result = make([][]float64, m.nao)
	for i := range result {
		result[i] = make([]float64, m.nao)
	}

	chan1 := make(chan CintData)
	maxGoroutines := runtime.GOMAXPROCS(-1)
	guard := make(chan struct{}, maxGoroutines)

	go func() {
		res_ptri := 0
		for i, _ := range m.Shells {
			res_ptrj := 0
			for j := 0; j < i+1; j++ {
				go func(ii, jj, prti, ptrj int) {
					chan1 <- m.CintKinetic(ii, jj, prti, ptrj)
					<-guard
				}(i, j, res_ptri, res_ptrj)
				res_ptrj += m.Shells[j]
			}
			res_ptri += m.Shells[i]
		}
	}()

	var wg sync.WaitGroup
	wg.Add(1)

	go func() {
		defer wg.Done()
		res_ptri := 0
		for i, _ := range m.Shells {
			res_ptrj := 0
			for j := 0; j < i+1; j++ {
				data := <-chan1

				c := 0
				for dj := 0; dj < data.dj; dj++ {
					for di := 0; di < data.di; di++ {
						result[data.ptr_i+di][data.ptr_j+dj] = float64(data.data[c])
						result[data.ptr_j+dj][data.ptr_i+di] = float64(data.data[c])
						c++
					}
				}
				res_ptrj += m.Shells[j]
			}
			res_ptri += m.Shells[i]
		}
	}()
	wg.Wait()

	fmt.Println("Kinetic energy integrals calculation done. Var 1.")
	return result
}

func (m *Molecule) ElecNuc() [][]float64 {
	var result [][]float64

	result = make([][]float64, m.nao)
	for i := range result {
		result[i] = make([]float64, m.nao)
	}

	chan1 := make(chan CintData)
	maxGoroutines := runtime.GOMAXPROCS(-1)
	guard := make(chan struct{}, maxGoroutines)

	go func() {
		res_ptri := 0
		for i, _ := range m.Shells {
			res_ptrj := 0
			for j := 0; j < i+1; j++ {
				go func(ii, jj, prti, ptrj int) {
					chan1 <- m.CintElecNuc(ii, jj, prti, ptrj)
					<-guard
				}(i, j, res_ptri, res_ptrj)
				res_ptrj += m.Shells[j]
			}
			res_ptri += m.Shells[i]
		}
	}()

	var wg sync.WaitGroup
	wg.Add(1)

	go func() {
		defer wg.Done()
		res_ptri := 0
		for i, _ := range m.Shells {
			res_ptrj := 0
			for j := 0; j < i+1; j++ {
				data := <-chan1

				c := 0
				for dj := 0; dj < data.dj; dj++ {
					for di := 0; di < data.di; di++ {
						result[data.ptr_i+di][data.ptr_j+dj] = float64(data.data[c])
						result[data.ptr_j+dj][data.ptr_i+di] = float64(data.data[c])
						c++
					}
				}
				res_ptrj += m.Shells[j]
			}
			res_ptri += m.Shells[i]
		}
	}()
	wg.Wait()

	fmt.Println("Electron-nuclei attraction integrals calculation done. Var 1.")
	return result
}

func (m *Molecule) CintElecElec(i, j, k, l, ptr_i, ptr_j, ptr_k, ptr_l int) Cint2Data {
	di := m.Shells[i]
	dj := m.Shells[j]
	dk := m.Shells[k]
	dl := m.Shells[l]
	buff_size := di * dj * dk * dl
	shls := [4]C.int{C.int(i), C.int(j), C.int(k), C.int(l)}
	buf := make([]C.double, buff_size)

	_ = C.cint2e_sph(&buf[0], &shls[0], &m.catm[0][0], m.natm, &m.cbas[0][0], m.nbas, &m.cenv[0], nil)

	return Cint2Data{di, dj, dk, dl, ptr_i, ptr_j, ptr_k, ptr_l, buf}
}

func (m *Molecule) ElecElec1() ([]int, []float64) {

	var res1 []int
	var res2 []float64

	chan1 := make(chan Cint2Data)
	maxGoroutines := runtime.GOMAXPROCS(-1)
	guard := make(chan struct{}, maxGoroutines)

	go func() {

		res_ptri := 0
		for i := 0; i < m.NShells; i++ {
			res_ptrj := 0
			for j := 0; j < i+1; j++ {
				res_ptrk := 0
				for k := 0; k < m.NShells; k++ {
					res_ptrl := 0
					for l := 0; l < k+1; l++ {
						go func(ii, jj, kk, ll, ptri, ptrj, ptrk, ptrl int) {
							chan1 <- m.CintElecElec(ii, jj, kk, ll, ptri, ptrj, ptrk, ptrl)
							<-guard
						}(i, j, k, l, res_ptri, res_ptrj, res_ptrk, res_ptrl)

						res_ptrl += m.Shells[l]
					}
					res_ptrk += m.Shells[k]
				}
				res_ptrj += m.Shells[j]
			}
			res_ptri += m.Shells[i]
		}
	}()

	var wg sync.WaitGroup
	wg.Add(1)

	go func() {
		defer wg.Done()
		count := 0
		res_ptri := 0
		for i := 0; i < m.NShells; i++ {
			res_ptrj := 0
			for j := 0; j < i+1; j++ {
				res_ptrk := 0
				for k := 0; k < m.NShells; k++ {
					res_ptrl := 0
					for l := 0; l < k+1; l++ {

						data := <-chan1
						c := 0
						for dl := 0; dl < data.dl; dl++ {
							for dk := 0; dk < data.dk; dk++ {
								for dj := 0; dj < data.dj; dj++ {
									for di := 0; di < data.di; di++ {
										if math.Abs(float64(data.data[c])-0) >= 1e-18 {
											if (data.ptr_j+dj <= data.ptr_i+di) && (data.ptr_l+dl <= data.ptr_k+dk){
											    res1 = append(res1, (data.ptr_i+di)*m.nao*m.nao*m.nao+
                                                                                                    (data.ptr_j+dj)*m.nao*m.nao+
                                                                                                    (data.ptr_k+dk)*m.nao+
                                                                                                    data.ptr_l+dl)
											    res2 = append(res2, float64(data.data[c]))
											}
											count += 1
										}
										c++
									}
								}
							}
						}
						res_ptrl += m.Shells[l]
					}
					res_ptrk += m.Shells[k]
				}
				res_ptrj += m.Shells[j]
			}
			res_ptri += m.Shells[i]
		}
		fmt.Println("\n*****Nonzero Vee integrals: ", count)
	}()
	wg.Wait()
	fmt.Println("Electron-Electron repulsion integrals calculation done. Var 2.")
	return res1, res2
}

func (m *Molecule) NucNuc() float64 {
	res := 0.0
	for i := range m.Atoms {
		for j := 0; j < i; j++ {
			res += float64(m.Atoms[i].Z) * float64(m.Atoms[j].Z) /
				math.Sqrt(math.Pow((m.Atoms[i].Coords[0]-m.Atoms[j].Coords[0])/a_B, 2)+
					math.Pow((m.Atoms[i].Coords[1]-m.Atoms[j].Coords[1])/a_B, 2)+
					math.Pow((m.Atoms[i].Coords[2]-m.Atoms[j].Coords[2])/a_B, 2))
		}
	}
	fmt.Println("Nuclei-Nuclei repulsion calculation done.")
	return res
}

func CalcContractionCoeffs(l int, c []float64) []float64 {
	res := make([]float64, len(c))
	for i := range res {
		res[i] = float64(C.CINTgto_norm(C.int(l), C.double(c[i])))
	}
	return res
}

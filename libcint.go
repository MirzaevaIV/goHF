// libcint.go --  This file is part of goHF project.
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

/*
#cgo CFLAGS: -g -Wall
#cgo LDFLAGS: ./cfiles/libcint.a
#include <stdlib.h>
#include <math.h>
#include "./cfiles/cint_funcs.h"
#include "./cfiles/cint.h"
#cgo LDFLAGS: -lm  -lquadmath
int cint1e_kin_sph(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env);
int cint1e_nuc_sph(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env);
int cint1e_ovlp_sph(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env);
//int cint2e_sph_optimizer(CINTOpt *opt, int *atm, int natm, int *bas, int nbas, double *env);
int cint2e_sph(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt);
*/
import "C"
import "math"

func (m *Molecule) Ovlp() [][]float64 {
	var result [][]float64
	cbas := CConvertInt(m.bas)
	catm := CConvertInt(m.atm)
	cenv := CConvertDouble(m.env)
	natm := C.int(len(m.atm))
	nbas := C.int(len(m.bas))
	var shls [2]C.int

	shells := make([]int, m.getNShells())
	for i := 0; i < len(shells); i++ {
		shells[i] = int(C.CINTcgto_spheric(C.int(i), &cbas[0][0]))
	}

	result = make([][]float64, m.nao)
	for i := range result {
		result[i] = make([]float64, m.nao)
	}

	res_ptri := 0

	for i, _ := range shells {
		res_ptrj := 0
		for j := 0; j < i+1; j++ {
			buff_size := int(shells[i] * shells[j])
			buf := make([]C.double, buff_size)
			shls[0] = C.int(i)
			shls[1] = C.int(j)

			_ = C.cint1e_ovlp_sph(&buf[0], &shls[0], &catm[0][0], natm, &cbas[0][0], nbas, &cenv[0])

			c := 0
			for dj := 0; dj < shells[j]; dj++ {
				for di := 0; di < shells[i]; di++ {
					result[res_ptri+di][res_ptrj+dj] = float64(buf[c])
					result[res_ptrj+dj][res_ptri+di] = float64(buf[c])
					c++
				}
			}
			res_ptrj += shells[j]
		}
		res_ptri += shells[i]
	}
	return result
}

func (m *Molecule) Kinetic() [][]float64 {
	var result [][]float64
	cbas := CConvertInt(m.bas)
	catm := CConvertInt(m.atm)
	cenv := CConvertDouble(m.env)
	natm := C.int(len(m.atm))
	nbas := C.int(len(m.bas))
	var shls [2]C.int

	shells := make([]int, m.getNShells())
	for i := 0; i < len(shells); i++ {
		shells[i] = int(C.CINTcgto_spheric(C.int(i), &cbas[0][0]))
	}

	result = make([][]float64, m.nao)
	for i := range result {
		result[i] = make([]float64, m.nao)
	}

	res_ptri := 0

	for i, _ := range shells {
		res_ptrj := 0
		for j := 0; j < i+1; j++ {
			buff_size := int(shells[i] * shells[j])
			buf := make([]C.double, buff_size)
			shls[0] = C.int(i)
			shls[1] = C.int(j)

			_ = C.cint1e_kin_sph(&buf[0], &shls[0], &catm[0][0], natm, &cbas[0][0], nbas, &cenv[0])

			c := 0
			for dj := 0; dj < shells[j]; dj++ {
				for di := 0; di < shells[i]; di++ {
					result[res_ptri+di][res_ptrj+dj] = float64(buf[c])
					result[res_ptrj+dj][res_ptri+di] = float64(buf[c])
					c++
				}
			}
			res_ptrj += shells[j]
		}
		res_ptri += shells[i]
	}
	return result
}

func (m *Molecule) ElecNuc() [][]float64 {
	var result [][]float64
	cbas := CConvertInt(m.bas)
	catm := CConvertInt(m.atm)
	cenv := CConvertDouble(m.env)
	natm := C.int(len(m.atm))
	nbas := C.int(len(m.bas))
	var shls [2]C.int

	shells := make([]int, m.getNShells())
	for i := 0; i < len(shells); i++ {
		shells[i] = int(C.CINTcgto_spheric(C.int(i), &cbas[0][0]))
	}

	result = make([][]float64, m.nao)
	for i := range result {
		result[i] = make([]float64, m.nao)
	}

	res_ptri := 0
	for i, _ := range shells {
		res_ptrj := 0
		for j := 0; j < i+1; j++ {
			buff_size := int(shells[i] * shells[j])
			buf := make([]C.double, buff_size)
			shls[0] = C.int(i)
			shls[1] = C.int(j)

			_ = C.cint1e_nuc_sph(&buf[0], &shls[0], &catm[0][0], natm, &cbas[0][0], nbas, &cenv[0])

			c := 0
			for dj := 0; dj < shells[j]; dj++ {
				for di := 0; di < shells[i]; di++ {
					result[res_ptri+di][res_ptrj+dj] = float64(buf[c])
					result[res_ptrj+dj][res_ptri+di] = float64(buf[c])

					c++
				}
			}
			res_ptrj += shells[j]
		}
		res_ptri += shells[i]
	}
	return result
}

func (m *Molecule) ElecElec() [][][][]float64 {
	var result [][][][]float64

	cbas := CConvertInt(m.bas)
	catm := CConvertInt(m.atm)
	cenv := CConvertDouble(m.env)
	natm := C.int(len(m.atm))
	nbas := C.int(len(m.bas))
	var shls [4]C.int

	nshells := m.getNShells()
	shells := make([]int, nshells)
	for i := 0; i < nshells; i++ {
		shells[i] = int(C.CINTcgto_spheric(C.int(i), &cbas[0][0]))
	}

	result = make([][][][]float64, m.nao)
	for i := range result {
		result[i] = make([][][]float64, m.nao)
		for j := range result[i] {
			result[i][j] = make([][]float64, m.nao)
			for k := range result[i][j] {
				result[i][j][k] = make([]float64, m.nao)
			}
		}
	}

	res_ptri := 0
	for i := 0; i < nshells; i++ {
		res_ptrj := 0
		for j := 0; j < i+1; j++ {
			res_ptrk := 0
			for k := 0; k < nshells; k++ {
				res_ptrl := 0
				for l := 0; l < k+1; l++ {
					shls[0] = C.int(i)
					shls[1] = C.int(j)
					shls[2] = C.int(k)
					shls[3] = C.int(l)

					buff_size := 1
					for sh_i := 0; sh_i < 4; sh_i++ {
						buff_size = buff_size * int(C.CINTcgto_spheric(shls[sh_i], &cbas[0][0]))
					}
					buf := make([]C.double, buff_size)

					_ = C.cint2e_sph(&buf[0], &shls[0], &catm[0][0], natm, &cbas[0][0], nbas, &cenv[0], nil)

					c := 0
					for dl := 0; dl < shells[l]; dl++ {
						for dk := 0; dk < shells[k]; dk++ {
							for dj := 0; dj < shells[j]; dj++ {
								for di := 0; di < shells[i]; di++ {
									result[res_ptri+di][res_ptrj+dj][res_ptrk+dk][res_ptrl+dl] = float64(buf[c])
									result[res_ptri+di][res_ptrj+dj][res_ptrl+dl][res_ptrk+dk] = float64(buf[c])
									result[res_ptrj+dj][res_ptri+di][res_ptrk+dk][res_ptrl+dl] = float64(buf[c])
									result[res_ptrj+dj][res_ptri+di][res_ptrl+dl][res_ptrk+dk] = float64(buf[c])
									c++
								}
							}
						}
					}
					res_ptrl += shells[l]
				}
				res_ptrk += shells[k]
			}
			res_ptrj += shells[j]
		}
		res_ptri += shells[i]
	}
	return result
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
	return res
}

func CalcContractionCoeffs(l int, c []float64) []float64 {
	res := make([]float64, len(c))
	for i := range res {
		res[i] = float64(C.CINTgto_norm(C.int(l), C.double(c[i])))
	}
	return res
}

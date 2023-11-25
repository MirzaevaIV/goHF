package main

// По мотивам https://github.com/nickelandcopper/HartreeFockPythonProgram/blob/main/Hartree_Fock_Program.ipynb
// По мотивам https://www.youtube.com/watch?v=eDAfpQIMde0&list=PL-hN8vfRaU7jSFHN1ZSAMNe_2nXhwAmzM

import (
	_ "encoding/csv"
	"fmt"
	_ "io"
	_ "log"
	"math"
	_ "os"
	_ "strings"
	_ "time"

	_ "gonum.org/v1/gonum"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/mathext"
)

type PrimitiveGaussian struct {
	Alpha  float64
	Coeff  float64
	Coords [3]float64 //Center coordinates
	L      [3]int     //Angular momentum vector
}

func (p PrimitiveGaussian) NormCoeff() float64 {
	return math.Pow((2 * p.Alpha / math.Pi), 0.75)
}

type AO struct {
	PGs []PrimitiveGaussian
}

type atom struct {
	Z      int
	Coords [3]float64
}

func QQ(v1, v2 [3]float64) float64 {
	vv1 := mat.NewVecDense(3, v1[:])
	vv2 := mat.NewVecDense(3, v2[:])
	dist := mat.NewVecDense(3, nil)
	dist.SubVec(vv2, vv1)
	dist.MulElemVec(dist, dist)
	res := mat.Sum(dist)
	return res
}

func CalcP(a1, a2 float64, v1, v2 [3]float64) [3]float64 {
	vv1 := mat.NewVecDense(3, v1[:])
	vv2 := mat.NewVecDense(3, v2[:])
	vres := mat.NewVecDense(3, nil)
	var res [3]float64
	vv1.ScaleVec(a1, vv1)
	vres.AddScaledVec(vv1, a2, vv2)
	for i, _ := range res {
		res[i] = vres.AtVec(i)
	}
	return res
}

func CalcPp(a float64, v [3]float64) [3]float64 {
	var res [3]float64
	for i, _ := range res {
		res[i] = v[i] / a
	}
	return res
}

func my_dot(v1, v2 [3]float64) float64 {
	vv1 := mat.NewVecDense(3, v1[:])
	vv2 := mat.NewVecDense(3, v2[:])
	res := mat.Dot(vv1, vv2)
	return res
}

func Overlap(m []AO) [][]float64 {
	N_AO := len(m)
	res := make([][]float64, N_AO)
	for i := range res {
		res[i] = make([]float64, N_AO)
	}

	for i := 0; i < N_AO; i++ {
		for j := 0; j < N_AO; j++ {
			for k := range m[i].PGs {
				for l := range m[j].PGs {
					N := m[i].PGs[k].NormCoeff() * m[j].PGs[l].NormCoeff()
					p := m[i].PGs[k].Alpha + m[j].PGs[l].Alpha
					q := m[i].PGs[k].Alpha * m[j].PGs[l].Alpha / p
					Q2 := QQ(m[i].PGs[k].Coords, m[j].PGs[l].Coords)

					res[i][j] += N * m[i].PGs[k].Coeff * m[j].PGs[l].Coeff * math.Exp(-q*Q2) * math.Pow((math.Pi/p), 1.5)

				}
			}
		}
	}
	return res

}

func Kinetic(m []AO) [][]float64 {
	N_AO := len(m)
	res := make([][]float64, N_AO)
	for i := range res {
		res[i] = make([]float64, N_AO)
	}

	for i := 0; i < N_AO; i++ {
		for j := 0; j < N_AO; j++ {
			for k := range m[i].PGs {
				for l := range m[j].PGs {
					N := m[i].PGs[k].NormCoeff() * m[j].PGs[l].NormCoeff()
					p := m[i].PGs[k].Alpha + m[j].PGs[l].Alpha
					q := m[i].PGs[k].Alpha * m[j].PGs[l].Alpha / p
					Q2 := QQ(m[i].PGs[k].Coords, m[j].PGs[l].Coords)

					P := CalcP(m[i].PGs[k].Alpha, m[j].PGs[l].Alpha, m[i].PGs[k].Coords, m[j].PGs[l].Coords)
					Pp := CalcPp(p, P)
					PGx2 := math.Pow((Pp[0] - m[j].PGs[l].Coords[0]), 2)
					PGy2 := math.Pow((Pp[1] - m[j].PGs[l].Coords[1]), 2)
					PGz2 := math.Pow((Pp[2] - m[j].PGs[l].Coords[2]), 2)

					c1c2 := m[i].PGs[k].Coeff * m[j].PGs[l].Coeff
					s := N * c1c2 * math.Exp(-q*Q2) * math.Pow((math.Pi/p), 1.5)

					res[i][j] += 3 * m[j].PGs[l].Alpha * s
					res[i][j] -= 2 * m[j].PGs[l].Alpha * m[j].PGs[l].Alpha * s * (PGx2 + 0.5/p)
					res[i][j] -= 2 * m[j].PGs[l].Alpha * m[j].PGs[l].Alpha * s * (PGy2 + 0.5/p)
					res[i][j] -= 2 * m[j].PGs[l].Alpha * m[j].PGs[l].Alpha * s * (PGz2 + 0.5/p)

				}
			}
		}
	}
	return res

}

func boys(x float64, n int) float64 {
	var res float64
	nf := float64(n)
	if x == 0 {
		res = 1.0 / (2.0*nf + 1)
	} else {
		res = mathext.GammaIncReg(nf+0.5, x) * math.Gamma(nf+0.5) * (1.0 / (2.0 * math.Pow(x, (nf+0.5))))
	}
	return res
}

func V_en(m []AO, a_list []atom) [][]float64 {
	N_AO := len(m)
	res := make([][]float64, N_AO)
	for i := range res {
		res[i] = make([]float64, N_AO)
	}

	for at := range a_list {
		for i := 0; i < N_AO; i++ {
			for j := 0; j < N_AO; j++ {
				for k := range m[i].PGs {
					for l := range m[j].PGs {
						N := m[i].PGs[k].NormCoeff() * m[j].PGs[l].NormCoeff()
						p := m[i].PGs[k].Alpha + m[j].PGs[l].Alpha
						q := m[i].PGs[k].Alpha * m[j].PGs[l].Alpha / p
						Q2 := QQ(m[i].PGs[k].Coords, m[j].PGs[l].Coords)

						P := CalcP(m[i].PGs[k].Alpha, m[j].PGs[l].Alpha, m[i].PGs[k].Coords, m[j].PGs[l].Coords)
						Pp := CalcPp(p, P)
						PGx2 := math.Pow((Pp[0] - a_list[at].Coords[0]), 2)
						PGy2 := math.Pow((Pp[1] - a_list[at].Coords[1]), 2)
						PGz2 := math.Pow((Pp[2] - a_list[at].Coords[2]), 2)
						PG2 := PGx2 + PGy2 + PGz2

						c1c2 := m[i].PGs[k].Coeff * m[j].PGs[l].Coeff

						res[i][j] += -float64(a_list[at].Z) * N * c1c2 * math.Exp(-q*Q2) * (2.0 * math.Pi / p) * boys(p*PG2, 0)

					}
				}
			}
		}
	}
	return res
}

func V_ee(m []AO) [][][][]float64 {
	N_AO := len(m)
	res := make([][][][]float64, N_AO) //create an empty array of required size
	for i := range res {
		res[i] = make([][][]float64, N_AO)
		for j := range res[i] {
			res[i][j] = make([][]float64, N_AO)
			for k := range res[i][j] {
				res[i][j][k] = make([]float64, N_AO)
			}
		}
	}

	for i := range m {
		for j := range m {
			for k := range m {
				for l := range m {
					for ii := range m[i].PGs {
						for jj := range m[j].PGs {
							for kk := range m[k].PGs {
								for ll := range m[l].PGs {
									N := m[i].PGs[ii].NormCoeff() * m[j].PGs[jj].NormCoeff() * m[k].PGs[kk].NormCoeff() * m[l].PGs[ll].NormCoeff()
									cicjckcl := m[i].PGs[ii].Coeff * m[j].PGs[jj].Coeff * m[k].PGs[kk].Coeff * m[l].PGs[ll].Coeff

									pij := m[i].PGs[ii].Alpha + m[j].PGs[jj].Alpha
									pkl := m[k].PGs[kk].Alpha + m[l].PGs[ll].Alpha
									Pij := CalcP(m[i].PGs[ii].Alpha, m[j].PGs[jj].Alpha, m[i].PGs[ii].Coords, m[j].PGs[jj].Coords)
									Pkl := CalcP(m[k].PGs[kk].Alpha, m[l].PGs[ll].Alpha, m[k].PGs[kk].Coords, m[l].PGs[ll].Coords)
									Ppij := CalcPp(pij, Pij)
									Ppkl := CalcPp(pkl, Pkl)
									PpijPpkl := [3]float64{(Ppij[0] - Ppkl[0]), (Ppij[1] - Ppkl[1]), (Ppij[2] - Ppkl[2])}
									PpijPpkl2 := math.Pow(PpijPpkl[0], 2) + math.Pow(PpijPpkl[1], 2) + math.Pow(PpijPpkl[2], 2)
									denom := (1.0 / pij) + (1.0 / pkl)

									qij := m[i].PGs[ii].Alpha * m[j].PGs[jj].Alpha / pij
									qkl := m[k].PGs[kk].Alpha * m[l].PGs[ll].Alpha / pkl
									Qij := [3]float64{(m[i].PGs[ii].Coords[0] - m[j].PGs[jj].Coords[0]), (m[i].PGs[ii].Coords[1] - m[j].PGs[jj].Coords[1]), (m[i].PGs[ii].Coords[2] - m[j].PGs[jj].Coords[2])}
									Qkl := [3]float64{(m[k].PGs[kk].Coords[0] - m[l].PGs[ll].Coords[0]), (m[k].PGs[kk].Coords[1] - m[l].PGs[ll].Coords[1]), (m[k].PGs[kk].Coords[2] - m[l].PGs[ll].Coords[2])}
									Q2ij := my_dot(Qij, Qij)
									Q2kl := my_dot(Qkl, Qkl)

									term1 := 2.0 * math.Pi * math.Pi / (pij * pkl)
									term2 := math.Sqrt(math.Pi / (pij + pkl))
									term3 := math.Exp(-qij * Q2ij)
									term4 := math.Exp(-qkl * Q2kl)

									res[i][j][k][l] += N * cicjckcl * term1 * term2 * term3 * term4 * boys(PpijPpkl2/denom, 0)

								}
							}
						}
					}
				}
			}
		}
	}
	return res
}

func E_nn(a_list []atom) float64 {
	res := 0.0
	for i := range a_list {
		for j := 0; j < i; j++ {
			res += float64(a_list[i].Z) * float64(a_list[j].Z) / math.Sqrt(math.Pow(a_list[i].Coords[0]-a_list[j].Coords[0], 2)+math.Pow(a_list[i].Coords[1]-a_list[j].Coords[1], 2)+math.Pow(a_list[i].Coords[2]-a_list[j].Coords[2], 2))
		}
	}
	return res
}

func SetSTO31G(dist float64) []AO {
	//STO-3G
	H1g1 := PrimitiveGaussian{0.3425250914e+01, 0.1543289673e+00, [3]float64{0.0, 0.0, 0.0}, [3]int{0, 0, 0}}
	H1g2 := PrimitiveGaussian{0.6239137298e+00, 0.5353281423e+00, [3]float64{0.0, 0.0, 0.0}, [3]int{0, 0, 0}}
	H1g3 := PrimitiveGaussian{0.1688554040e+00, 0.4446345422e+00, [3]float64{0.0, 0.0, 0.0}, [3]int{0, 0, 0}}

	H2g1 := PrimitiveGaussian{0.3425250914e+01, 0.1543289673e+00, [3]float64{dist, 0.0, 0.0}, [3]int{0, 0, 0}}
	H2g2 := PrimitiveGaussian{0.6239137298e+00, 0.5353281423e+00, [3]float64{dist, 0.0, 0.0}, [3]int{0, 0, 0}}
	H2g3 := PrimitiveGaussian{0.1688554040e+00, 0.4446345422e+00, [3]float64{dist, 0.0, 0.0}, [3]int{0, 0, 0}}

	H1 := AO{[]PrimitiveGaussian{H1g1, H1g2, H1g3}}
	H2 := AO{[]PrimitiveGaussian{H2g1, H2g2, H2g3}}

	res := []AO{H1, H2}
	return res
}

func Set631G(dist float64) []AO {
	//6-31+G
	H1g1 := PrimitiveGaussian{0.1873113696e+02, 0.3349460434e-01, [3]float64{0.0, 0.0, 0.0}, [3]int{0, 0, 0}}
	H1g2 := PrimitiveGaussian{0.2825394365e+01, 0.2347269535e+00, [3]float64{0.0, 0.0, 0.0}, [3]int{0, 0, 0}}
	H1g3 := PrimitiveGaussian{0.6401216923e+00, 0.8137573261e+00, [3]float64{0.0, 0.0, 0.0}, [3]int{0, 0, 0}}
	H1g4 := PrimitiveGaussian{0.1612777588e+00, 1.0000000, [3]float64{0.0, 0.0, 0.0}, [3]int{0, 0, 0}}

	H2g1 := PrimitiveGaussian{0.1873113696e+02, 0.3349460434e-01, [3]float64{dist, 0.0, 0.0}, [3]int{0, 0, 0}}
	H2g2 := PrimitiveGaussian{0.2825394365e+01, 0.2347269535e+00, [3]float64{dist, 0.0, 0.0}, [3]int{0, 0, 0}}
	H2g3 := PrimitiveGaussian{0.6401216923e+00, 0.8137573261e+00, [3]float64{dist, 0.0, 0.0}, [3]int{0, 0, 0}}
	H2g4 := PrimitiveGaussian{0.1612777588e+00, 1.0000000, [3]float64{dist, 0.0, 0.0}, [3]int{0, 0, 0}}

	H1_1s := AO{[]PrimitiveGaussian{H1g1, H1g2, H1g3}}
	H1_2s := AO{[]PrimitiveGaussian{H1g4}}
	H2_1s := AO{[]PrimitiveGaussian{H2g1, H2g2, H2g3}}
	H2_2s := AO{[]PrimitiveGaussian{H2g4}}

	res := []AO{H1_1s, H1_2s, H2_1s, H2_2s}
	return res
}

func SetAtoms(dist float64) []atom {
	H1_at := atom{1, [3]float64{0.0, 0.0, 0.0}}
	H2_at := atom{1, [3]float64{dist, 0.0, 0.0}}
	res := []atom{H1_at, H2_at}
	return res
}

type MolParameters struct {
	atomlist []atom
	Nocc     int
	mol      []AO
	S        [][]float64
	T        [][]float64
	Ven      [][]float64
	Vee      [][][][]float64
}

func compute_G(den_mat [][]float64, Vee [][][][]float64) [][]float64 {
	res := make([][]float64, len(den_mat))
	for i := range res {
		res[i] = make([]float64, len(den_mat))
	}
	for i := range Vee {
		for j := range Vee {
			for k := range Vee {
				for l := range Vee {
					J := Vee[i][j][k][l]
					K := Vee[i][l][k][j]
					res[i][j] += den_mat[k][l] * (J - 0.5*K)
				}
			}
		}
	}
	//fmt.Println("debug mode: denmat, vee, G")
	//fmt.Println(den_mat)
	//fmt.Println(Vee)
	//fmt.Println(res)
	return res
}

func compute_density_mat(MOs [][]float64, N_occ int) [][]float64 {
	n_basis := len(MOs)
	res := make([][]float64, n_basis)
	for i := range res {
		res[i] = make([]float64, n_basis)
	}
	occ := 2.0
	for i := 0; i < n_basis; i++ {
		for j := 0; j < n_basis; j++ {
			for oo := 0; oo < N_occ; oo++ {
				C := MOs[i][oo]
				C_dag := MOs[j][oo]
				res[i][j] += occ * C * C_dag
			}
		}
	}
	return res
}

func compute_energy(D, T, Vne, G [][]float64) float64 {
	//fmt.Println("D, T, Vne, G")
	//fmt.Println(D)
	//fmt.Println(T)
	//fmt.Println(Vne)
	//fmt.Println(G)
	n_basis := len(T)
	res := 0.0
	for i := 0; i < n_basis; i++ {
		for j := 0; j < n_basis; j++ {
			res += D[i][j] * (T[i][j] + Vne[i][j] + 0.5*G[i][j])
			//fmt.Println(i, j, res)
		}
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

func scf_cycle(MP MolParameters, Tol float64, MaxSteps int) float64 {
	res := 0.0
	E_prev := 0.0
	n_basis := len(MP.mol)
	density_mat := make([][]float64, n_basis)
	for i := range density_mat {
		density_mat[i] = make([]float64, n_basis)
	}
	H1 := mat.NewDense(n_basis, n_basis, flatten(MP.T)) //one-electron term
	H1.Add(H1, mat.NewDense(n_basis, n_basis, flatten(MP.Ven)))
	// 1. Enter into the SCF cycles
	for i := 0; i < MaxSteps; i++ {
		E_prev = res

		// 2. Compute the 2 electron term
		G := compute_G(density_mat, MP.Vee)

		// 3. Form initial F
		F := mat.NewDense(n_basis, n_basis, flatten(G))
		F.Add(F, H1)

		//Diagonalize and normalize right side of equation.  S^{-1/2} S S^{-1/2}
		Smat := mat.NewSymDense(n_basis, flatten(MP.S))
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

		//Transform the left side. S^{-1/2} F S^{-1/2}
		F.Mul(F, &SSqrtInv)
		F.Mul(&SSqrtInv, F)
		FSym := mat.NewSymDense(n_basis, F.RawMatrix().Data)
		ok = eigsym.Factorize(FSym, true)
		if !ok {
			fmt.Println("Transformed F eigendecomposition failed")
		}

		eigsym.VectorsTo(&ev)
		ev.Mul(&SSqrtInv, &ev)

		//make new MO coefficients matrix
		MOs := make([][]float64, n_basis)
		for i := range MOs {
			MOs[i] = ev.RawRowView(i)
		}

		// 4. Form new density matrix using MOs
		density_mat = compute_density_mat(MOs, MP.Nocc)

		// 5. Compute electronic_energy expectation value (res)
		res = compute_energy(density_mat, MP.T, MP.Ven, G)
		// 6. Check convergence
		if math.Abs(E_prev-res) < Tol {
			fmt.Println("SCF converged after step ", i)
			//fmt.Println(density_mat)
			return res
		}
	}
	fmt.Println("Warning! SCF NOT converged after step ", MaxSteps)
	return res
}

func main() {

	dist := 0.6
	results := [][]float64{}
	step := 0.1
	NSteps := 20
	Tolerance := 1e-8
	SCFMaxSteps := 20
	N_occ := 1
	for i := 0; i <= NSteps; i++ {
		atomlist := SetAtoms(dist)
		mol := Set631G(dist)
		S := Overlap(mol)
		T := Kinetic(mol)
		Ven := V_en(mol, atomlist)
		Vee := V_ee(mol)
		MolPars := MolParameters{atomlist, N_occ, mol, S, T, Ven, Vee}
		E := scf_cycle(MolPars, Tolerance, SCFMaxSteps)
		fmt.Println(MolPars.S, MolPars.T, MolPars.Ven, MolPars.Vee)
		Enn := E_nn(atomlist)
		results = append(results, []float64{dist, Enn, E + Enn})
		dist += step
	}
	for i := range results {
		fmt.Println(results[i])
	}
}

// molecule.go --  This file is part of goHF project.
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

import (
	"strconv"
	"strings"

	"golang.org/x/exp/slices"
)

type Molecule struct {
	Atoms      []Atom
	atm        [][]int
	bas        [][]int
	env        []float64
	nao, Nelec int
}

func (m *Molecule) addAtoms(data []string, start int, end int) {
	for i := start; i < end+1; i++ {
		var atm Atom
		words := strings.Fields(data[i])
		atm.Z = slices.Index(ElemData.Symb, words[0])
		atm.Name = words[0] + strconv.Itoa(1+i-start)
		if len(words) > 3 {
			x, _ := strconv.ParseFloat(words[1], 64)
			y, _ := strconv.ParseFloat(words[2], 64)
			z, _ := strconv.ParseFloat(words[3], 64)
			atm.Coords = [3]float64{x, y, z}
		} else {
			ErrorLogger.Println("Incorrect format of coordinates for atom: ", atm.Name)
		}
		m.Atoms = append(m.Atoms, atm)
	}
}

func (m *Molecule) getBasis(bName string) {
	bFile := "./data/basis/" + strings.ToLower(strings.Fields(bName)[0]) + ".txt"
	data, err := ReadFileLines(bFile)
	if err != nil {
		ErrorLogger.Println("Cannot read basis file ", bFile, ": ", err)
	}
	for i, atm := range m.Atoms {
		for j, str := range data {
			words := strings.Fields(str)
			if len(words) > 1 {
				if (len(words[0])>2) && (words[1] == strings.ToUpper(ElemData.Symb[atm.Z])) {
					OutputLogger.Println(i+1, "Basis for atom ", atm.Name, ": ", data[j+1])
					m.Atoms[i].getBasis(data, j+2)
				}
			}
		}
	}
}

func (m *Molecule) getNelec() int {
	result := 0
	for _, a := range m.Atoms {
		result += a.Z
	}
	return result
}

func (m *Molecule) getNShells() int {
	result := 0
	for _, a := range m.Atoms {
		result += len(a.Basis)
	}
	return result
}

func (m *Molecule) buildCINTdata() {
	m.atm = [][]int{}
	m.bas = [][]int{}
	m.env = make([]float64, 20)
	ptr_env := 20
	for _, a := range m.Atoms {
		m.atm = append(m.atm, []int{a.Z, ptr_env, 0, 0, 0, 0})
		m.env = append(m.env, a.Coords[0]/a_B, a.Coords[1]/a_B, a.Coords[2]/a_B)
		ptr_env += 3
	}
	for i, a := range m.Atoms {
		for _, o := range a.Basis {
			var exps []float64
			var preExps []float64
			for _, pg := range o.Funcs {
				exps = append(exps, pg.zeta)
				preExps = append(preExps, pg.preExp)
			}
			contr_coef := CalcContractionCoeffs(o.l, exps)
			for _, el := range exps {
				m.env = append(m.env, el)
			}
			for idx, el := range preExps {
				m.env = append(m.env, el*contr_coef[idx])
			}
			m.bas = append(m.bas, []int{i, o.l, o.nPrim, 1, 0, ptr_env, ptr_env + o.nPrim, 0})
			ptr_env += 2 * o.nPrim
			m.nao += (2*o.l + 1)
		}
	}
}

func (mol *Molecule) CalculateIntegrals() ([][]float64, [][]float64, [][]float64, [][][][]float64) {
	mol.buildCINTdata()
	//OutputLogger.Println(mol)
	//printOutputDelimiter()

	//printOutputDelimiter()
	//OutputLogger.Println("Overlap Integrals: ")
	//printOutputDelimiter()
	S := mol.Ovlp()
	//OutputLogger.Println(S)

	//printOutputDelimiter()
	//OutputLogger.Println("Kinetic Energy Integrals: ")
	//printOutputDelimiter()
	T := mol.Kinetic()
	//OutputLogger.Println(T)

	//printOutputDelimiter()
	//OutputLogger.Println("Nuclei Potential Energy Integrals: ")
	//printOutputDelimiter()
	Vn := mol.ElecNuc()
	//OutputLogger.Println(Vn)

	//printOutputDelimiter()
	//OutputLogger.Println("2-electron integrals: ")
	Vee := mol.ElecElec()
	//OutputLogger.Print(len(Vee))

	//printOutputDelimiter()

	return S, T, Vn, Vee
}

func (atm *Atom) getBasis(data []string, pos int) {
	nOrbs, _ := strconv.Atoi(strings.Fields(data[pos])[0])
	pos++
	for k := 0; k < nOrbs; k++ {
		var orb Orbital
		orb.n, _ = strconv.Atoi(strings.Fields(data[pos])[0])
		orb.l, _ = strconv.Atoi(strings.Fields(data[pos])[1])
		orb.nPrim, _ = strconv.Atoi(strings.Fields(data[pos])[2])
		pos++
		for l := 0; l < orb.nPrim; l++ {
			var pg PrimitiveGauss
			pg.zeta, _ = strconv.ParseFloat(strings.Fields(data[pos])[0], 64)
			pg.preExp, _ = strconv.ParseFloat(strings.Fields(data[pos])[1], 64)
			orb.Funcs = append(orb.Funcs, pg)
			pos++
		}
		atm.Basis = append(atm.Basis, orb)
	}
}

type Atom struct {
	Z      int
	Name   string
	Coords [3]float64
	Basis  []Orbital
}

type Orbital struct {
	n, l, nPrim int
	Funcs       []PrimitiveGauss
}

type PrimitiveGauss struct {
	zeta, preExp float64
}

type Mendeleev struct {
	Z          []int
	Symb, Name []string
	Mass       []float64
}

func (m *Mendeleev) build() {
	data, err := ReadFileLines("./data/mendeleev.csv")
	if err != nil {
		ErrorLogger.Println("Cannot read elements database file './data/mendeleev.csv' file: ", err)
	}
	for i, str := range data {
		if i == 0 {
		} else {
			words := strings.Split(str, ",")
			z, _ := strconv.Atoi(words[0])
			mass, _ := strconv.ParseFloat(words[3], 64)
			m.Z = append(m.Z, z)
			m.Mass = append(m.Mass, mass)
			m.Symb = append(m.Symb, words[1])
			m.Name = append(m.Name, words[2])
		}
	}
}

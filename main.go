// main.go --  This file is part of goHF project.
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
	"log"
	_ "math"
	"os"
	"runtime"
	"strings"
	"strconv"
	_ "time"

	_ "golang.org/x/exp/slices"
	_ "gonum.org/v1/gonum/mat"
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
				OutputLogger.Print("Parsing input. Atoms block found at lines ", atom_start, " -- ", atom_end, ".")
			}
			if strings.ToLower(words[0]) == "basis" {
				basis = true
				basisName = data[i+1]
				_ = findBlockEnd(i, data, "Basis")
				OutputLogger.Print("Parsing input. Basis block found.", basisName)
			}
			if strings.ToLower(words[0]) == "nprocs" {
				nprocs, _ := strconv.Atoi(words[1])
                runtime.GOMAXPROCS(nprocs)
				OutputLogger.Print("Parsing input. Number of threads set to "+words[1]+".")
			}
		}
	}
	if !atoms {
		ErrorLogger.Fatal("Parsing input. No Atoms found.")
	} else {
		mol.addAtoms(data, atom_start+1, atom_end-1)
	}
	if !basis {
		OutputLogger.Println("Parsing input. No Basis found. Using default basis: STO-3G.")
	} else {
		mol.getBasis(basisName)
		mol.setNPrims()
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


func main() {
	runtime.GOMAXPROCS(1)

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

	mol := processInput(inpData) //--- create molecule from input data

	calc := mol.RHFinit() // -- initialize RHF structure
	
	calc.BuildDensMat()
	
    EEE := calc.SCF_DIIS()

	OutputLogger.Println("Nuclei Repulsion Energy: ", mol.NucNuc(), " a.u.")
	printOutputDelimiter()
	fmt.Println("Nuc energy = ", mol.NucNuc(), " a.u.")

	fmt.Println("Final total energy = ", EEE, " a.u.")
	OutputLogger.Println("Final total energy = ", EEE, " a.u.")
	printOutputDelimiter()

	MyMemDebug()

	OutputLogger.Println("\n")
	InfoLogger.Println("Exiting goHF...")
	fmt.Println("goHF done.")
}

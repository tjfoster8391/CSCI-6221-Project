// RecursiveNoise.go
package main

import (
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/png"
	"math/rand"
	"os"

	//"time"
	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/app"
	"fyne.io/fyne/v2/canvas"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/widget"
)

var termList map[int]PerlinTerm // map to list perlin terms
var mapI int
var termsPrinted int
var length int // make this a power of two
var imageArray [512][512]float64

var originalNoise [512][512]float64 // store original to use for zooming in to start new recursion

var pixelThreshold int

var intiated bool
var zoomedIn bool
var imagePath string

var myImage image.RGBA

var w fyne.Window
var content *fyne.Container

var thisImage *canvas.Image

//var thisColor color.Color

// make initial
type PerlinTerm struct {
	period    int // period, 1/f of that term
	amplitude float64
	colors    [512][512]float64 // interpolated colors
	lattice   [512][512]float64 // lattice of random values

	// length int?
	// for zooming in and out
}

type coord struct {
	x float64
	y float64
}

func main() {

	// create GUI components

	a := app.New()

	w = a.NewWindow("Recursive Noise")

	zoomIn := widget.NewButton("zoom In", func() { zoomIn() }) // zooms in
	initB := widget.NewButton("init", func() { initF() })      // calls initialization, creating the first term and recursively creating more

	content = container.NewBorder(zoomIn, initB, nil, nil, nil)

	w.SetContent(content)

	w.Resize(fyne.NewSize(512, 552)) // approximte dimensions, with extra for buttons

	w.ShowAndRun()

	/*fmt.Println("zoomIn called")
	fmt.Println()

	zoomIn()*/
	//
}

func initF() {

	//fmt.Println("initF")
	length = 512
	zoomedIn = false

	imagePath = "Noise.png"

	pixelThreshold = 1.0
	termList = make(map[int]PerlinTerm) // create map to store terms, with ordered keys using "mapI"
	firstTerm := createTerm(*new(PerlinTerm), 32, 0.5)
	mapI = 0
	termList[mapI] = *firstTerm
	mapI++
	addTermsRecursive(*firstTerm) // creates all the terms until the period is less than the size of a pixel
	sumTerms()

	drawNoise()
	intiated = true
}
func drawNoise() { // draws array containing summation of perlin terms
	fmt.Println("noise drawn")

	myImage := image.NewRGBA(image.Rect(0, 0, 512, 512)) // create drawing object to draw noise
	//thisColor := color.RGBA{0, 0, 0, 255}
	myFile, err := os.Create(imagePath) // creates output image
	if err != nil {
		panic(err)
	}
	backUpFile, err := os.Create("OriginalNoise.png") // backup image to compare to zoomed in one
	if err != nil {
		panic(err)
	}

	/*fmt.Println("noise being drawn")
	fmt.Printf("\n")*/
	for i := 0; i < len(imageArray); i++ {
		for j := 0; j < len(imageArray[i]); j++ {
			greyCoeff := imageArray[i][j] // getting greyscale value from image array

			colorValue := int(greyCoeff * 255.0) // convert it to the 0 to 255 scale

			thisColor := color.RGBA{uint8(colorValue), uint8(colorValue), uint8(colorValue), 255} // create a color for this point
			currentDot := image.Rect(i, j, i+1, j+1)                                              // must create 2x2 pixels cause the draw function requires a rectangle

			draw.Draw(myImage, currentDot, &image.Uniform{thisColor}, image.ZP, draw.Src)

		}
		//fmt.Printf("\n")
	}
	//fmt.Println("image updated: ", termsPrinted)
	png.Encode(myFile, myImage) // create file
	if !zoomedIn {
		png.Encode(backUpFile, myImage) // create additional file
	}
	defer myFile.Close()

	thisImage = canvas.NewImageFromFile(imagePath) // create image that can be put in container

	thisImage.ScaleMode = canvas.ImageScaleSmooth

	content.Add(thisImage) // add it

}

func interpolate(x int, y int, x1 int, x2 int, y1 int, y2 int, fx1y1 float64, fx1y2 float64, fx2y1 float64, fx2y2 float64) float64 {
	// takes nearest values from the lattice and bilinear interpolation to them.
	fxy1 := (float64(x2-x)/float64(x2-x1))*fx1y1 + (float64(x-x1)/float64(x2-x1))*fx2y1

	fxy2 := (float64(x2-x)/float64(x2-x1))*fx1y2 + (float64(x-x1)/float64(x2-x1))*fx1y2

	fxy := (float64(y2-y)/float64(y2-y1))*fxy1 + (float64(y-y1)/float64(y2-y1))*fxy2

	///fmt.Println("interpolation return: ", fxy)
	return fxy

}

func addTermsRecursive(term PerlinTerm) {
	//fmt.Println("add terms recursive called ")

	if (term.period / 2) < pixelThreshold { // if a term's period is smaller than a pixel, no need to keep making more
		//fmt.Println("base case reached")
		var newTerm PerlinTerm = *createTerm(term, term.period, term.amplitude)
		termList[mapI] = newTerm
		mapI++
		return
	} else {
		var newTerm PerlinTerm = *createTerm(term, term.period/2, term.amplitude/2) // term created, added to map, recusion continued
		termList[mapI] = newTerm
		mapI++
		addTermsRecursive(newTerm)
	}
}

func createTerm(previousTerm PerlinTerm, period int, amplitude float64) *PerlinTerm {
	// must rely on previous term
	// all the previous
	//fmt.Println("create term called ")
	var newTerm PerlinTerm // new term and variables created for this term.
	var lattice = [512][512]float64{}
	var colors = [512][512]float64{}
	newTerm.lattice = lattice
	newTerm.colors = colors
	newTerm.period = period
	newTerm.amplitude = amplitude

	//fmt.Println("new term: period: ", period, " ", "amplitude", amplitude)

	// do I need this in to
	if period == 32 { // creating the first term
		for i := 0; i < length; i += period {
			for j := 0; j < length; j += period { // lattice points generated randominly
				newTerm.lattice[i][j] = rand.Float64()
				if newTerm.lattice[i][j] > 0.5 { // make sure in range from 0 to 0.5
					newTerm.lattice[i][j] /= 2
				}
				newTerm.colors[i][j] = newTerm.lattice[i][j] // color gets the same points to reference the lattice
				originalNoise[i][j] = newTerm.lattice[i][j]

			}
		}

	} else {
		for i := 0; i < length-1; i += period { // if not first term, shrink previous term array and copy it into the 4 quadrants of the new array
			// self similarity
			for j := 0; j < length-1; j += period {
				//fmt.Println("index check: ", i, j)
				newTerm.colors[i/2][j/2] = previousTerm.lattice[i][j] / 2
				newTerm.lattice[i/2][j/2] = previousTerm.lattice[i][j] / 2

				newTerm.colors[i/2+length/2][j/2] = previousTerm.lattice[i][j] / 2
				newTerm.lattice[i/2+length/2][j/2] = previousTerm.lattice[i][j] / 2

				newTerm.colors[i/2][j/2+length/2] = previousTerm.lattice[i][j] / 2
				newTerm.lattice[i/2][j/2+length/2] = previousTerm.lattice[i][j] / 2

				newTerm.colors[i/2+length/2][j/2+length/2] = previousTerm.lattice[i][j] / 2
				newTerm.lattice[i/2+length/2][j/2+length/2] = previousTerm.lattice[i][j] / 2
				// the halved period version of the previous lattice takes up one quarter of the new one

			}
		}
	}
	//fmt.Println("period: ", newTerm.period)
	for k := 0; k < len(newTerm.colors); k += 2 { // interpolate in color array based on lattice points regardless of term
		for l := 0; l < len(newTerm.colors[k]); l += 2 {
			if newTerm.lattice[k][l] == 0 {
				// lattice coordinates surrounding the point, to be interpolated

				prevX := k - (k % newTerm.period) // calulate surrounding values
				nextX := k + ((k + newTerm.period - ((k + newTerm.period) % newTerm.period)) - k)
				prevY := l - (l % newTerm.period)
				nextY := l + ((l + newTerm.period - ((l + newTerm.period) % newTerm.period)) - l)

				if nextY > 511 { // prevent index out of bounds errors.
					nextY = 511
				}
				if nextX > 511 {
					nextX = 511
				}

				newTerm.colors[k][l] = interpolate(k, l, prevX, nextX, nextY, prevY, newTerm.lattice[prevX][nextY],
					newTerm.lattice[prevX][prevY], newTerm.lattice[nextX][nextY], newTerm.lattice[nextX][prevY])

			}

		}
	}

	return &newTerm // return time to add to map

}

func zoomIn() {
	if intiated == false {
		fmt.Println("initialize first")
	} else {
		//fmt.Println("zoom in called")
		zoomedIn = true // other methods have a difference
		mapI = 0        // restart map
		t := termList[0]

		t.period *= 2

		var replaceLattice [512][512]float64 // put top left corner of original as entirety of new one.
		var replaceColors [512][512]float64
		for i := 0; i < len(t.lattice); i += t.period {
			for j := 0; j < len(t.lattice); j += t.period {
				replaceLattice[i][j] = originalNoise[i/2][j/2]
				replaceColors[i][j] = originalNoise[i/2][j/2]
			}
		}
		t.colors = replaceColors
		t.lattice = replaceLattice

		for k := 0; k < len(t.colors); k += 2 { // interpolate
			for l := 0; l < len(t.colors[k]); l += 2 {
				if t.lattice[k][l] == 0 {

					prevX := k - (k % t.period)
					nextX := k + ((k + t.period - ((k + t.period) % t.period)) - k)
					prevY := l - (l % t.period)
					nextY := l + ((l + t.period - ((l + t.period) % t.period)) - l)

					if nextY > 511 {
						nextY = 511
					}
					if nextX > 511 {
						nextX = 511
					} // avoid index out of bounds

					t.colors[k][l] = interpolate(k, l, prevX, nextX, nextY, prevY, t.lattice[prevX][nextY],
						t.lattice[prevX][prevY], t.lattice[nextX][nextY], t.lattice[nextX][prevY])

				}

			}
		}

		termList[mapI] = t
		mapI++

		addTermsRecursive(termList[0])

		sumTerms()
		drawNoise()
	}
}

func sumTerms() {

	var total [512][512]float64
	//fmt.Println("terms summed")
	//fmt.Println()

	//fmt.Println("mapI: ", mapI)
	termsPrinted = 1 // track terms in original generation

	//fmt.Println("mapI: ", mapI, "length of list", len(termList))
	for k := 0; k <= len(termList); k++ {
		//fmt.Println("k: ", k, ". termList element", termList[k].period, "first point: ")
		if !zoomedIn {
			termsPrinted++
		}
		for i := 0; i < length; i++ {
			for j := 0; j < length; j++ {

				total[i][j] += termList[k].colors[i][j]
				//fmt.Print(termList[k].colors[i][j], " ")
				imageArray = total // might be using pointers wrong

			}
			//fmt.Println()
		}
		if !zoomedIn {
			drawNoise()
			//fmt.Println("new term printed: ", termsPrinted)
		} else {
			if !(k < termsPrinted) { // if zoomed in, draw amount of terms from original all at once, draw new ones incrementally
				drawNoise()
				//fmt.Println("new term printed: ", termsPrinted)
			}
		}
		//fmt.Println()
	}

}

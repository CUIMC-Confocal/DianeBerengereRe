//Creating circular well area, cut it into tiles and select only in-focus tiles with normalized forcus socre
//with adaptive cell detection based on DAPI channel: mean
//with adaptive threshold to classify Green+Red+ cells
//Yanyan Chen (09/18/2026; Columbia University; yc4569@cumc.columbia.edu)

//Import OpenCV and Qupath OpenCV library functions for image processing

import org.bytedeco.opencv.global.opencv_imgproc
import org.bytedeco.opencv.global.opencv_core
import org.bytedeco.opencv.opencv_core.Mat
import org.bytedeco.opencv.opencv_core.Size
import qupath.opencv.tools.OpenCVTools
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
import qupath.lib.objects.PathObjects

//Define Functions (4 functions here)

//Function: FocusScore calculation: using the variance of Laplacian

def calculateFocusScoreForTiles() {

    def imageData = getCurrentImageData()
    def server = imageData.getServer()
    def tiles = getAnnotationObjects()

    for (tile in tiles) {

        def roi = tile.getROI()
        def request = RegionRequest.createInstance(server.getPath(), 1.0, roi)
        def img = server.readRegion(request)

        def mat = OpenCVTools.imageToMat(img)

        // Channels: green channel
        Mat green = new Mat()
        //Mat blue = new Mat()

        //opencv_core.extractChannel(mat, blue, 0)
        opencv_core.extractChannel(mat, green, 1)

        // Combine blue and green channels
        //Mat combined = new Mat()
        //opencv_core.addWeighted(
        //   green, 0.5,
        //   blue, 0.5,
        //   0,
        //   combined
        //)

        // Convert green channel to floating point
        Mat greenFloat = new Mat()
        green.convertTo(
            greenFloat,
            opencv_core.CV_64F
        )

        // Calculate mean intensity
        Mat intensityMean = new Mat()
        Mat intensityStd = new Mat()

        opencv_core.meanStdDev(
            greenFloat,
            intensityMean,
            intensityStd
        )

        double meanIntensity =
            intensityMean.createBuffer().get(0)

        // Blur
        opencv_imgproc.GaussianBlur(
            greenFloat,
            greenFloat,
            new Size(3,3),
            0
        )

        // Laplacian
        Mat lap = new Mat()

        opencv_imgproc.Laplacian(
            greenFloat,
            lap,
            opencv_core.CV_64F
        )

        // Calculate Laplacian variance
        Mat lapMean = new Mat()
        Mat lapStd = new Mat()

        opencv_core.meanStdDev(
            lap,
            lapMean,
            lapStd
        )

        double lapStdValue =
            lapStd.createBuffer().get(0)

        double lapVariance =
            lapStdValue * lapStdValue

        // Normalize focus score by mean intensity x1000000 does not change the ranking but easier to see with a larger number
        double normalizedFocusScore = 0

        if (meanIntensity > 0) {
            normalizedFocusScore =
                (lapVariance / (meanIntensity * meanIntensity))*1000000
        }

        // Save measurements
        tile.getMeasurementList().put(
            "MeanIntensity_Green",
            meanIntensity
        )

        tile.getMeasurementList().put(
            "FocusScore_Original_Green",
            lapVariance
        )

        tile.getMeasurementList().put(
            "FocusScore_Normalized_Green",
            normalizedFocusScore
        )
    }

    fireHierarchyUpdate()

    println "FocusScore calculated for " +
            tiles.size() + " tiles"
}

//Function: Select in-focus tiles that has focus score higher than a threshold value

def markInFocusTiles(double thresholdValue) {
    def tiles = getAnnotationObjects()

    // Create the InFocus class
    def goodClass = getPathClass("InFocus") ?: addPathClass("InFocus")

    def inFocusTiles = []

    tiles.each { tile ->
        def score = tile.getMeasurementList().get("FocusScore_Normalized_Green")
        if (score != null && score > thresholdValue) {
            tile.setPathClass(goodClass)
            inFocusTiles << tile
        }
    }

    fireHierarchyUpdate()

    return inFocusTiles
}

//Function: calculate mean and std of the cell signals (Green/Red channel) from the detected cells in the inFocus tiles, will be used in the next function

def calculateStatistics(objects, measurementName) {

    def values = objects.collect { obj ->
        obj.getMeasurementList().get(measurementName)
    }.findAll { value ->
        value != null &&
        !Double.isNaN(value) &&
        !Double.isInfinite(value)
    }

    if (values.isEmpty()) {
        return null
    }

    double mean = values.sum() / values.size()

    double sumSquared = values.collect {
        (it - mean) * (it - mean)
    }.sum()

    double sd = values.size() > 1 ?
        Math.sqrt(sumSquared / (values.size() - 1)) :
        0.0

    return [
        values: values,
        mean: mean,
        sd: sd
    ]
}

// Function: Classify Green/Red positive cells,  using only cells detected in the in-focus tiles
// Statistics are calculated separately for each image

def classifyGreenRedCells(inFocusTiles) {

    // Get all cell detections from the in-focus tiles only

    def inFocusDetections = inFocusTiles.collectMany { tile ->

        tile.getChildObjects().findAll { obj ->
            obj.isDetection()
        }

    }.unique()


    // Check that detections exist

    if (inFocusDetections.isEmpty()) {
        println "No cell detections found inside inFocusTiles!"
        return
    }


    println "Cell detections used for classification: " +
            inFocusDetections.size()


    // Green statistics

    def greenMeasurement = "Cell: Channel 2 mean"

    def greenStats = calculateStatistics(
        inFocusDetections,
        greenMeasurement
    )

    if (greenStats == null) {
        println "No valid green measurements found!"
        return
    }

    double greenMean = greenStats.mean
    double greenSD   = greenStats.sd

    // 1+ = mean
    // 2+ = mean + 1 sd
    // 3+ = mean + 2 sd

    double greenThreshold1 = greenMean
    double greenThreshold2 = greenMean + greenSD
    double greenThreshold3 = greenMean + 2 * greenSD


    println "Green signal intensity"
    println "Mean = " + greenMean
    println "SD   = " + greenSD
    println "1+ threshold = " + greenThreshold1
    println "2+ threshold = " + greenThreshold2
    println "3+ threshold = " + greenThreshold3


    // Red statistics

    def redMeasurement = "Cell: Channel 1 mean"

    def redStats = calculateStatistics(
        inFocusDetections,
        redMeasurement
    )

    if (redStats == null) {
        println "No valid red measurements found!"
        return
    }

    double redMean = redStats.mean
    double redSD   = redStats.sd

    // 1+ = mean
    // 2+ = mean + 1 sd
    // 3+ = mean + 2 sd

    double redThreshold1 = redMean
    double redThreshold2 = redMean + redSD
    double redThreshold3 = redMean + 2 * redSD


    println "Red signal intensity"
    println "Mean = " + redMean
    println "SD   = " + redSD
    println "1+ threshold = " + redThreshold1
    println "2+ threshold = " + redThreshold2
    println "3+ threshold = " + redThreshold3


    // Create classes

    def negativeClass = getPathClass("Negative")
    def greenClass    = getPathClass("Green+")
    def redClass      = getPathClass("Red+")
    def greenRedClass = getPathClass("Green+Red+")


    // Classify cells

    int greenCount     = 0
    int redCount       = 0
    int greenRedCount  = 0
    int greenOnlyCount = 0
    int redOnlyCount   = 0
    int negativeCount  = 0


    inFocusDetections.each { cell ->

        def greenValue =
            cell.getMeasurementList().get(greenMeasurement)

        def redValue =
            cell.getMeasurementList().get(redMeasurement)

        // Skip cells with missing measurements

        if (greenValue == null ||
            redValue == null ||
            Double.isNaN(greenValue) ||
            Double.isNaN(redValue)) {

            return
        }


        // Determine Green and Red positive cells
        // Currently using threshold 1+ = image-specific mean 

        boolean isGreen = greenValue >= greenThreshold1
        boolean isRed   = redValue >= redThreshold1


        // Determine Green+ and Red+

        if (isGreen && isRed) {

            cell.setPathClass(greenRedClass)

            greenRedCount++
            greenCount++
            redCount++
        }


        // Determine Green+ only

        else if (isGreen) {

            cell.setPathClass(greenClass)

            greenOnlyCount++
            greenCount++
        }


        // Determine Red+ only

        else if (isRed) {

            cell.setPathClass(redClass)

            redOnlyCount++
            redCount++
        }


        // Determine the negative cells (nerither Green+ or Red+)

        else {

            cell.setPathClass(negativeClass)

            negativeCount++
        }
    }

    // Update QUPath

    fireHierarchyUpdate()
}

//-----------------------------------------------------------------
//Main batch pipeline

def project = getProject()

//Entry in project.getImageList().take(6) means only run the first 6 images
for (entry in project.getImageList()) {

    //Print the file name that is going to be processed
    println "Processing: " + entry.getImageName()   

    def imageData = entry.readImageData()
    setBatchProjectAndImage(project, imageData)

    //Reset and select the image type
    clearAllObjects()
    setImageType('FLUORESCENCE')
    resetSelection()
    
    // Set pixel size in um
    double pixelsize= 1.0714
    
    // Set center coordinates (x, y) and size (diameter in pixels)
    double x = 3952 / pixelsize
    double y = 3860 / pixelsize
    double diameter = 6500 /pixelsize
    
    // Top-left corner for the bounding box of the circle
    double rx = x - (diameter / 2)
    double ry = y - (diameter / 2)
    
    // Get default image plane (z=0, t=0)
    def plane = ImagePlane.getDefaultPlane()
    
    // Create ellipse ROI with equal width and height for a circle
    def roi = ROIs.createEllipseROI(rx, ry, diameter, diameter, plane)
    
    // Create and add the annotation object to the current hierarchy
    def annotation = PathObjects.createAnnotationObject(roi)

    imageData.getHierarchy().addObject(annotation)

    // Get Channel 3 (DAPI) mean and standard deviation
    selectAnnotations();
    runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons":2.0,"region":"ROI","tileSizeMicrons":25.0,"channel1":true,"channel2":true,"channel3":true,"doMean":true,"doStdDev":true,"doMinMax":true,"doMedian":true,"doHaralick":false,"haralickMin":NaN,"haralickMax":NaN,"haralickDistance":1,"haralickBins":32}')


    // Get Channel 3 mean and standard deviation
    def measurements = annotation.getMeasurementList()

    def channel3Mean = measurements.get("ROI: 2.00 µm per pixel: Channel 3: Mean")
    def channel3Std = measurements.get("ROI: 2.00 µm per pixel: Channel 3: Std.dev.")

    // Threshold = mean + standard deviation, channel3Mean + channel3Std
    // Threshold = mean
    def threshold_CellDetction = channel3Mean     

    println "Channel 3 Mean = " + channel3Mean
    println "Channel 3 Std = " + channel3Std
    println "Watershed Threshold = " + threshold_CellDetction


   //Cut the WellArea into tiles and remove the parent annotation
    selectObjects(getAnnotationObjects())
    runPlugin('qupath.lib.algorithms.TilerPlugin',
        '{"tileSizeMicrons":1450.0,"trimToROI":true,"makeAnnotations":true,"removeParentAnnotation":true}')

    entry.saveImageData(imageData)


    //Calculate the FocusScore for each tile
    calculateFocusScoreForTiles()

    // Filter in-focus tiles based on the threshold. Note: threshold = 4000 is an empirical number 
    double threshold = 4000
    def inFocusTiles = markInFocusTiles(threshold)

    println "In-focus tiles: " + inFocusTiles.size()

    if (inFocusTiles.isEmpty()) {
        println "No in-focus tiles → skipping image"
        continue
    }


    // Cell detection only on the InFocus tiles
    selectObjects(inFocusTiles)

    //Do the cell detection based on the DAPI channel, the threshold is adpat to each image: it's Mean
    runPlugin(
    'qupath.imagej.detect.cells.WatershedCellDetection',
    '{"detectionImage":"Channel 3","requestedPixelSizeMicrons":0.5,"backgroundRadiusMicrons":8.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":2.0,"minAreaMicrons":20.0,"maxAreaMicrons":250.0,"threshold":' + threshold_CellDetction + ',"watershedPostProcess":true,"cellExpansionMicrons":4.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true}'
    )

    //Found out the green neuron cells, cells with red signal and green neuron with red signal 

    classifyGreenRedCells(inFocusTiles)

    // Save results to the results folder under the Project, each image will have a ImageName.cvs file 
    def resultsDir = buildFilePath(PROJECT_BASE_DIR, "results")
    mkdirs(resultsDir)

    def outPath = buildFilePath(resultsDir, entry.getImageName() + ".csv")

    saveAnnotationMeasurements(outPath)

    entry.saveImageData(imageData)

}

println "Done!"
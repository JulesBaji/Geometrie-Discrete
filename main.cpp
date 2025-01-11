#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageSelector.h>
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/ColorBrightnessColorMap.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include "DGtal/io/Color.h"

#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/Color.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/helpers/StdDefs.h"

#include <DGtal/geometry/curves/GreedySegmentation.h>
#include <DGtal/geometry/curves/ArithmeticalDSSComputer.h>
#include <DGtal/geometry/curves/FreemanChain.h>
#include <sys/stat.h>

#include <vector>
#include <string>
#include <iostream>
#include <filesystem>

using namespace std;
using namespace DGtal;
using namespace Z2i;

template<class T>
Curve getBoundary(T & object, std::vector<Point> boundaryPoints)
{
    //Khalimsky space
    KSpace kSpace;
    // we need to add a margine to prevent situations such that an object touch the bourder of the domain
    kSpace.init( object.domain().lowerBound() - Point(1,1),
                   object.domain().upperBound() + Point(1,1), true );

    // 1) Call Surfaces::findABel() to find a cell which belongs to the border
    Z2i::SCell bel = Surfaces<Z2i::KSpace>::findABel(kSpace, object.pointSet(), 10000);
    // 2) Call Surfece::track2DBoundaryPoints to extract the boundary of the object
    SurfelAdjacency<KSpace::dimension> SAdj( true );
    Surfaces<KSpace>::track2DBoundaryPoints( boundaryPoints, kSpace, SAdj, object.pointSet(), bel );
    Curve boundaryCurve;
    // 3) Create a curve from a vector
    boundaryCurve.initFromVector(boundaryPoints);
    return boundaryCurve;
}

template<class T>
void sendToBoard( Board2D & board, T & p_Object, DGtal::Color p_Color) {
    board << CustomStyle( p_Object.className(), new DGtal::CustomFillColor(p_Color));
    board << p_Object;
}

template<class T>
void polygonizeBoundary(T & boundary, std::vector<LibBoard::Point> & libboardPoints)
{
    typedef FreemanChain<int> Contour4;
    typedef ArithmeticalDSSComputer<Contour4::ConstIterator, int, 4> DSS4;
    typedef GreedySegmentation<DSS4> Decomposition4;

    Z2i::Curve::PointsRange range = boundary.getPointsRange();
    std::string fmc; 
    Z2i::Point prevPoint = *range.begin(); 
    Z2i::Point beginPoint = prevPoint;
    Z2i::Point currentPoint;
    std::vector<Z2i::Point> directions = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};

    fmc += std::to_string(beginPoint[0]) + " " + std::to_string(beginPoint[1]) + " ";
    
    for (auto it = std::next(range.begin()); it != range.end(); ++it) {
        currentPoint = *it;
        Z2i::Point delta = currentPoint - prevPoint;

        for (int i = 0; i < directions.size(); ++i) {
            if (delta == directions[i]) {
                fmc += std::to_string(i);
                break;
            }
        }

        prevPoint = currentPoint;
    }

    std::stringstream ss(stringstream::in | stringstream::out);
    ss << fmc << std::endl;
    Contour4 theContour(ss);

    Decomposition4 theDecomposition(theContour.begin(), theContour.end(), DSS4());
    LibBoard::Point closureVertex(0, 0);
    for (Decomposition4::SegmentComputerIterator it = theDecomposition.begin(), itEnd = theDecomposition.end(); it != itEnd; ++it){
        if (closureVertex.x == 0) {
            closureVertex.x = (*it).Uf()[0];
            closureVertex.y = (*it).Uf()[1];
        }
        libboardPoints.push_back(LibBoard::Point((*it).Uf()[0], (*it).Uf()[1]));
    }   
    libboardPoints.push_back(LibBoard::Point(closureVertex.x, closureVertex.y));
}

double calculatePolygonArea(const std::vector<LibBoard::Point>& polygon) {
    int n = polygon.size();
    if (n < 3) return 0.0;

    double area = 0.0;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n; 
        area += polygon[i].x * polygon[j].y - polygon[i].y * polygon[j].x;
    }
    return std::abs(area) / 2.0;
}

double calculatePolygonPerimeter(const std::vector<LibBoard::Point>& polygon) {
    double perimeter = 0.0;
    int n = polygon.size();
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n; // Indice du point suivant
        int dx = polygon[j].x - polygon[i].x;
        int dy = polygon[j].y - polygon[i].y;
        perimeter += std::sqrt(dx * dx + dy * dy);
    }
    return perimeter;
}

double calculateMean(const std::vector<double>& values) {
    double sum = 0;
    for (double v : values) {
        sum += v;
    }
    return sum / values.size();
}

double calculateVariance(const std::vector<double>& values, double mean) {
    double variance = 0;
    for (double v : values) {
        variance += (v - mean) * (v - mean);
    }
    return variance / values.size();
}

void calculateAndDisplayStatistics(
    const std::string& title,
    const std::vector<double>& values)
{
    if (values.empty())
    {
        std::cout << title << ": No data available.\n";
        return;
    }

    double mean = calculateMean(values);
    double variance = calculateVariance(values, mean);
    std::cout << title << ":" << std::endl;
    std::cout << "Mean = " << mean << ", Variance = " << variance << std::endl;
}

int main(int argc, char** argv)
{
    // Step 1
    setlocale(LC_NUMERIC, "us_US"); // To prevent French local settings
    typedef ImageSelector<Domain, unsigned char>::Type Image; // type of image
    typedef DigitalSetSelector<Domain, BIG_DS + HIGH_BEL_DS>::Type DigitalSet; // Digital set type
    typedef Object<DT8_4, DigitalSet> ObjectType; // Digital object type

    // List of image paths
    std::vector<std::string> imagePaths = {
        "../RiceGrains/Rice_camargue_seg_bin.pgm",
        "../RiceGrains/Rice_japonais_seg_bin.pgm",
        "../RiceGrains/Rice_basmati_seg_bin.pgm"
    };

    std::vector<double> camargueAreas, camarguePolygonAreas, camarguePerimeters, camarguePolygonPerimeters, camargueCircularities;
    std::vector<double> japonaisAreas, japonaisPolygonAreas, japonaisPerimeters, japonaisPolygonPerimeters, japonaisCircularities;
    std::vector<double> basmatiAreas, basmatiPolygonAreas, basmatiPerimeters, basmatiPolygonPerimeters, basmatiCircularities;

    // Iterate over each image
    for (size_t imgIdx = 0; imgIdx < imagePaths.size(); ++imgIdx)
    {
        const std::string& imagePath = imagePaths[imgIdx];
        // Read image
        Image image = PGMReader<Image>::importPGM(imagePath.c_str()); // provide image path
        Z2i::DigitalSet set2d(image.domain());

        // Populate the digital set from the image
        SetFromImage<Z2i::DigitalSet>::append<Image>(set2d, image, 1, 255);

        // Create a vector to store objects (connected components)
        std::vector<ObjectType> objects;
        std::back_insert_iterator<std::vector<ObjectType>> inserter(objects);
        
        // Set the adjacency pair and obtain the connected components
        ObjectType object8_4(dt8_4, set2d);
        object8_4.writeComponents(inserter);
        std::cout << "Image: " << imagePath << " - Number of components: " << objects.size() << std::endl;

        // Step 2: Removal of partially outside components
        objects.erase(
            std::remove_if(objects.begin(), objects.end(), [&](const ObjectType& obj) {
                const auto& domain = image.domain();
                const auto& lowerBound = domain.lowerBound();
                const auto& upperBound = domain.upperBound();

                return std::any_of(obj.pointSet().begin(), obj.pointSet().end(), [&](const Z2i::Point& p) {
                    return !(p[0] > lowerBound[0] && p[0] < upperBound[0] &&
                             p[1] > lowerBound[1] && p[1] < upperBound[1]);
                });
            }),
            objects.end()
        );
        std::cout << "Image: " << imagePath << " - Number of entirely inside components: " << objects.size() << std::endl;

        // Create image-specific folders
        std::string imageName = imagePath.substr(imagePath.find_last_of("/") + 1, imagePath.find_last_of(".") - imagePath.find_last_of("/") - 1);
        std::string boundaryDir = "boundary_" + imageName;
        std::string componentDir = "component_" + imageName;
        std::string polygonizationDir = "polygonization_" + imageName;

        mkdir(boundaryDir.c_str(), 0777);
        mkdir(componentDir.c_str(), 0777);
        mkdir(polygonizationDir.c_str(), 0777);

        // Iterate through each object and process them
        for (size_t i = 0; i < objects.size(); ++i)
        {
            // Step 2
            Board2D aBoard;
            sendToBoard(aBoard, objects[i], Color::Red);
            std::string componentFilename = componentDir + "/out_component_" + std::to_string(i) + ".pdf";
            aBoard.saveCairo(componentFilename.c_str(), Board2D::CairoPDF);

            // Step 3
            Board2D aBoard2;
            std::vector<Point> boundaryPoints;
            Curve boundaryCurve = getBoundary(objects[i], boundaryPoints);
            sendToBoard(aBoard2, boundaryCurve, Color::Red);
            std::string boundaryFilename = boundaryDir + "/out_boundary_" + std::to_string(i) + ".pdf";
            aBoard2.saveCairo(boundaryFilename.c_str(), Board2D::CairoPDF);

            // Step 4
            std::vector<LibBoard::Point> libboardPoints;
            polygonizeBoundary(boundaryCurve, libboardPoints);
            Board2D aBoard3;
            aBoard3.drawPolyline(libboardPoints);
            std::string polygonizationFilename = polygonizationDir + "/polygonization_" + std::to_string(i) + ".pdf";
            aBoard3.saveCairo(polygonizationFilename.c_str(), Board2D::CairoPDF);

            // Step 5
            double area = objects[i].pointSet().size();
            double polygonArea = calculatePolygonArea(libboardPoints);

            // Step 6
            double perimeter = boundaryCurve.size() ;
            double polygonPerimeter = calculatePolygonPerimeter(libboardPoints);

            // Step 7
            double circularity = (4 * M_PI * polygonArea) / (polygonPerimeter * polygonPerimeter);

            if (imgIdx == 0) { // Rice Camargue
                camargueAreas.push_back(area);
                camarguePolygonAreas.push_back(polygonArea);
                camarguePerimeters.push_back(perimeter);
                camarguePolygonPerimeters.push_back(polygonPerimeter);
                camargueCircularities.push_back(circularity);
            } else if (imgIdx == 1) { // Rice Japonais
                japonaisAreas.push_back(area);
                japonaisPolygonAreas.push_back(polygonArea);
                japonaisPerimeters.push_back(perimeter);
                japonaisPolygonPerimeters.push_back(polygonPerimeter);
                japonaisCircularities.push_back(circularity);
            } else if (imgIdx == 2) { // Rice Basmati
                basmatiAreas.push_back(area);
                basmatiPolygonAreas.push_back(polygonArea);
                basmatiPerimeters.push_back(perimeter);
                basmatiPolygonPerimeters.push_back(polygonPerimeter);
                basmatiCircularities.push_back(circularity);
            }
        }
    }

    std::cout << "Results for Rice Camargue:" << std::endl;
    calculateAndDisplayStatistics("Average and Variance for Areas", camargueAreas);
    calculateAndDisplayStatistics("Average and Variance for Polygonal Areas", camarguePolygonAreas);
    calculateAndDisplayStatistics("Average and Variance for Perimeters", camarguePerimeters);
    calculateAndDisplayStatistics("Average and Variance for Polygonal Perimeters", camarguePolygonPerimeters);
    calculateAndDisplayStatistics("Average and Variance for Circularity", camargueCircularities);

    std::cout << "Results for Rice Japonais:" << std::endl;
    calculateAndDisplayStatistics("Average and Variance for Areas", japonaisAreas);
    calculateAndDisplayStatistics("Average and Variance for Polygonal Areas", japonaisPolygonAreas);
    calculateAndDisplayStatistics("Average and Variance for Perimeters", japonaisPerimeters);
    calculateAndDisplayStatistics("Average and Variance for Polygonal Perimeters", japonaisPolygonPerimeters);
    calculateAndDisplayStatistics("Average and Variance for Circularity", japonaisCircularities);

    std::cout << "Results for Rice Basmati:" << std::endl;
    calculateAndDisplayStatistics("Average and Variance for Areas", basmatiAreas);
    calculateAndDisplayStatistics("Average and Variance for Polygonal Areas", basmatiPolygonAreas);
    calculateAndDisplayStatistics("Average and Variance for Perimeters", basmatiPerimeters);
    calculateAndDisplayStatistics("Average and Variance for Polygonal Perimeters", basmatiPolygonPerimeters);
    calculateAndDisplayStatistics("Average and Variance for Circularity", basmatiCircularities);

    return 0;
}
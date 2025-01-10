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

using namespace std;
using namespace DGtal;
using namespace Z2i;

template<class T>
std::vector<Z2i::Point> getBoundary(T & object)
{
    //Khalimsky space
    KSpace kSpace;
    // we need to add a margine to prevent situations such that an object touch the bourder of the domain
    kSpace.init( object.domain().lowerBound() - Point(1,1),
                   object.domain().upperBound() + Point(1,1), true );

    // 1) Call Surfaces::findABel() to find a cell which belongs to the border
    Z2i::SCell cell_bel = Surfaces<Z2i::KSpace>::findABel(kSpace, object.pointSet(), 10000);
    // 2) Call Surfece::track2DBoundaryPoints to extract the boundary of the object
    Surfaces<Z2i::KSpace>::track2DBoundaryPoints(kSpace, object.pointSet(), cell_bel, boundaryPoints);
    return boundaryPoints;
}

template<class T>
void sendToBoard( Board2D & board, T & p_Object, DGtal::Color p_Color) {
    board << CustomStyle( p_Object.className(), new DGtal::CustomFillColor(p_Color));
    board << p_Object;
}

template<class T>
void polygonizeBoundary(T & boundary, std::vector<Z2i::Point> & polylinePoints)
{
    typedef FreemanChain<int> Contour4;
    typedef ArithmeticalDSSComputer<Contour4::ConstIterator, int, 4> DSS4;
    typedef GreedySegmentation<DSS4> Decomposition4;

    // Convertir la frontière en une chaîne de Freeman
    std::stringstream ss;
    if (boundary.size() > 0) {
        // Ajouter le premier point de la frontière
        ss << boundary[0][0] << " " << boundary[0][1] << " ";  // Premier point

        // Générer les déplacements élémentaires pour chaque point suivant
        for (size_t i = 1; i < boundary.size(); ++i) {
            int dx = boundary[i][0] - boundary[i - 1][0];
            int dy = boundary[i][1] - boundary[i - 1][1];

            // Ajouter un code de direction basé sur les déplacements
            if (dx == 1 && dy == 0) ss << "1";  // Droite
            else if (dx == -1 && dy == 0) ss << "3";  // Gauche
            else if (dx == 0 && dy == 1) ss << "2";  // Haut
            else if (dx == 0 && dy == -1) ss << "0";  // Bas
            else if (dx == 1 && dy == 1) ss << "4";  // Haut-Droite
            else if (dx == -1 && dy == 1) ss << "5";  // Haut-Gauche
            else if (dx == 1 && dy == -1) ss << "7";  // Bas-Droite
            else if (dx == -1 && dy == -1) ss << "6";  // Bas-Gauche
        }
    }

    // Afficher la chaîne de Freeman pour vérification
    std::cout << "Freeman chain: " << ss.str() << std::endl;

    // Créer la chaîne de Freeman
    Contour4 theContour(ss);

    // Appliquer la segmentation greedy pour détecter les segments droits numériques
    Decomposition4 theDecomposition(theContour.begin(), theContour.end(), DSS4());

    // Ajouter les points de la segmentation à polylinePoints
    for (auto it = theDecomposition.begin(); it != theDecomposition.end(); ++it) {
        // Ajouter les points de chaque segment à polylinePoints
        for (auto pt : it->primitive().points()) {
            polylinePoints.push_back(pt);
        }
    }
}

double calculatePolygonArea(const std::vector<Z2i::Point>& polygon) {
    int n = vertices.size();
    if (n < 3) return 0.0;

    double area = 0.0;
    int n = polygon.size();
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n; 
        area += polygon[i][0] * polygon[j][1] - polygon[i][1] * polygon[j][0];
    }
    return std::abs(area) / 2.0;
}

double calculatePolygonPerimeter(const std::vector<Z2i::Point>& polygon) {
    double perimeter = 0.0;
    int n = polygon.size();
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;  // Le point suivant, avec retour au début du tableau
        // Calcul de la distance entre deux points consécutifs
        int dx = polygon[j][0] - polygon[i][0];
        int dy = polygon[j][1] - polygon[i][1];
        perimeter += std::sqrt(dx * dx + dy * dy);
    }
    return perimeter;
}

int main(int argc, char** argv)
{
    // Step 1
    setlocale(LC_NUMERIC, "us_US"); //To prevent French local settings
    typedef ImageSelector<Domain, unsigned char >::Type Image; // type of image
    typedef DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSet; // Digital set type
    typedef Object<DT8_4, DigitalSet> ObjectType; // Digital object type
    // typedef Object<DT4_8, DigitalSet> ObjectType; // Digital object type

    // read an image
    Image image = PGMReader<Image>::importPGM ("../RiceGrains/Rice_japonais_seg_bin.pgm"); // you have to provide a correct path as a parameter
    // 1) make a "digital set" of proper size
    Z2i::DigitalSet set2d (image.domain());
    // 2) populate a digital set from the image using SetFromImage::append()
    SetFromImage<Z2i::DigitalSet>::append<Image>(set2d, image, 1, 255);
    // 3) Create a digital object from the digital set
    std::vector< ObjectType > objects;          // All connected components are going to be stored in it
    std::back_insert_iterator< std::vector< ObjectType > > inserter( objects ); // Iterator used to populated "objects".
    // 4) Set the adjacency pair and obtain the connected components using "writeComponents"
    // ObjectType object4_8(dt4_8, set2d);
    ObjectType object8_4(dt8_4, set2d);
    object8_4.writeComponents(inserter);
    std::cout << " number of components : " << objects.size() << endl; // Right now size of "objects" is the number of conected components

    // Step 2
    objects.erase(
        std::remove_if(objects.begin(), objects.end(), [&](const ObjectType& obj) {
            const auto& domain = image.domain();
            const auto& lowerBound = domain.lowerBound();
            const auto& upperBound = domain.upperBound();

            // Vérifie si tous les points de l'objet sont dans les limites
            return std::any_of(obj.pointSet().begin(), obj.pointSet().end(), [&](const Z2i::Point& p) {
                return !(p[0] > lowerBound[0] && p[0] < upperBound[0] &&
                        p[1] > lowerBound[1] && p[1] < upperBound[1]);
            });
        }),
        objects.end()
    );
    std::cout << "number of entire components = " << objects.size() << std::endl;

    Board2D aBoard;
    sendToBoard(aBoard, objects[0], Color::Red);   // send the connected component "objects[0]" to "aBoard"
    aBoard.saveCairo("out.pdf",Board2D::CairoPDF); // do not forget to change the path!
    
    // Step 3
    Board2D aBoard2;
    std::vector<Z2i::Point> boundaryPoints = getBoundary(objects[0]);
    Curve boundaryCurve;
    // Create a curve from a vector
    for (const auto& point : boundaryPoints) {
        //curve.initFromSCellsVector(boundryVector.at(i));
        boundaryCurve.push_back(point);
    }
    sendToBoard(aBoard2, boundaryCurve, Color::Red);
    aBoard2.saveCairo("out_boundary_" + std::to_string(0) + ".pdf", Board2D::CairoPDF);

    // Step 4
    std::vector<Z2i::Point> polylinePoints;
    polygonizeBoundary(boundaryCurve, polylinePoints);
    // Visualiser la polygonisation
    Board2D aBoard3;
    aBoard3.drawPolyline(polylinePoints);
    aBoard3.saveCairo("poygonisation_" + std::to_string(0) + ".pdf", Board2D::CairoPDF);

    // Step 5
    std::cout << "Aire grille pour l'objet " << 0 << " : " << objects[0].pointSet().size() << std::endl;
    double polygonArea = calculatePolygonArea(polylinePoints);
    std::cout << "Aire polygone pour l'objet " << 0 << " : " << polygonArea << std::endl;

    // Step 6
    std::cout << "Périmètre (1-cellules) pour l'objet " << 0 << " : " << boundaryCurve.size(); << std::endl;
    double polygonPerimeter = calculatePolygonPerimeter(polylinePoints);
    std::cout << "Périmètre du polygone pour l'objet " << 0 << " : " << polygonPerimeter << std::endl;

    // Step 7
    double circularity  = (4 * M_PI * polygonArea) / (polygonPerimeter * polygonPerimeter);
    std::cout << "Circularité pour l'objet " << 0 << " : " << circularity << std::endl;

    return 0;
}
#include <iostream>
#include <vector>
#include <cstdlib> // per rand() e srand()
#include <ctime> // per time()
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

#define THRASHOLD_SEGENMENTATION_DERIVATIVE 50
#define WINDOW_SIZE_SMOOTH 1

struct Point {
    float x, y;     // coordinates

    Point(float r, float th) : 
        // conversion from polar to cartesian
        x(r*cos(th)), // x = radiuscos(theta)
        y(r*sin(th)) // y = radiussin(theta)
    {}
    
    Point(float x_input, float y_input,bool is_cartesian) : 
        // conversion from polar to cartesian
        x(x_input), // x = radiuscos(theta)
        y(y_input) // y = radiussin(theta)
    {}
};

std::vector<std::vector<Point>> circle_centers(std::vector<float> ranges, float incrementTheta,float thetaBegin);
std::vector<float> clean_first_20(std::vector<float> ranges);
void divide_in_segment(std::vector<float> ranges,std::vector<float> f_deriv, std::vector<std::vector<float>> &segmentedData, std::vector<int> &segment_index_begin);
void calc_first_derivative_polar(std::vector<float> data, std::vector<float> &output, float incrementTheta);
void segmentPolar_to_segmentCartesian(std::vector<std::vector<float>> segmentedData, std::vector<std::vector<Point>> &segmentedDataCartesian, float initial_theta, float incrementTheta);
std::vector<std::vector<Point>>smooth_data_cartesian(std::vector<std::vector<Point>> dataSegemntedPoint, int window_size);
std::vector<std::vector<Point>>calc_first_derivative_segment_smooth(std::vector<std::vector<Point>> data);


float findiff(float fxh, float fx, float h) {
    return (fxh - fx) / h;
}
int main() {
    //inizio debug lettura file
    //inizializza il generatore di numeri casuali
    std::vector<float> distanze;
    //leggi il file di nome data.txt
    std::string nome_file = "data.txt";

    // Apri il file in modalità di lettura
    std::ifstream file_input(nome_file);

    if (!file_input.is_open()) {
        std::cerr << "Impossibile aprire il file: " << nome_file << std::endl;
        return 1; // Termina il programma con un codice di errore
    }

    // Leggi il contenuto del file
    std::string linea;
    while (std::getline(file_input, linea)) {
        // Utilizza un stringstream per analizzare i numeri da ogni riga
        std::stringstream stream_linea(linea);
        
        std::string numero_str;
        while (std::getline(stream_linea, numero_str, ',')) {
            try {
                // Converte la stringa in un numero double e lo aggiunge al vettore
                double numero = std::stof(numero_str);
                distanze.push_back(numero);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Errore nella conversione del numero: o"<< numero_str <<"o " << e.what()<< std::endl;
            }
        }
    }

        // Chiudi il file dopo aver finito di leggere
        file_input.close();
    //fine debug lettura file
    //stampo le distanze
    //calolca i centri dei cerchi
    float incrementTheta = 0.005774, thetaBegin = -1.9198600053787231;
    std::cout << "incremento theta: " << incrementTheta << "\n";
    std::vector<std::vector<Point>> circle_centers_var = circle_centers(distanze, incrementTheta,thetaBegin);
    std::cout << "size: " << circle_centers_var.size() << "\n";
    //stampa i cerchi nel terminale
    for(int i = 0; i < circle_centers_var.size(); i++){
        std::cout << "segmento: " << i << "\n";
        for(int j = 0; j < circle_centers_var[i].size(); j++){
            std::cout << "x: " << circle_centers_var[i][j].x << " y: " << circle_centers_var[i][j].y << "\n";
        }
    }
   
    return 0;
}










std::vector<float> clean_first_20(std::vector<float> ranges){
    std::vector<float> ranges_clean;
    for(int i = 20; i < ranges.size()-20; i++){
        ranges_clean.push_back(ranges[i]);
    }
    return ranges_clean;
}
void calc_first_derivative_polar(std::vector<float> data, std::vector<float> &output, float incrementTheta){
    for(int i = 0; i < data.size()-1; i++){
        output.push_back(findiff(data[i+1], data[i], incrementTheta));
    }
}

void divide_in_segment(std::vector<float> ranges,std::vector<float> f_deriv, std::vector<std::vector<float>> &segmentedData, std::vector<int> &segment_index_begin){
    int segment_index = 0;
    double flag_begin = true;

    std::vector<float> segment;
    for(int i = 0; i < ranges.size()-1; i++){
        if(std::fabs(f_deriv[i]) < THRASHOLD_SEGENMENTATION_DERIVATIVE){
            if(flag_begin){
                segment_index_begin.push_back(i);
                flag_begin = false;
            }
            segment.push_back(ranges[i]);
        }
        else{
            if(segment.size() > 0)
                segmentedData.push_back(segment);
            segment.clear();
            flag_begin = true;
        }
    }
    if(segment.size() > 0)
        segmentedData.push_back(segment);
}

void segmentPolar_to_segmentCartesian(std::vector<std::vector<float>> segmentedData, std::vector<std::vector<Point>> &segmentedDataCartesian, float initial_theta, float incrementTheta){
    for(int i = 0; i < segmentedData.size(); i++){
        std::vector<Point> segmentCartesian;
        float initial_theta_i = initial_theta + (segmentedData[i].size()/2)*incrementTheta;
        std::cout << "initial_theta_i: " << initial_theta_i << "\n";
        for(int j = 0; j < segmentedData[i].size(); j++){
            segmentCartesian.push_back(Point(segmentedData[i][j], initial_theta_i + (j*incrementTheta)));
        }
        segmentedDataCartesian.push_back(segmentCartesian);
    }
}

std::vector<std::vector<Point>>smooth_data_cartesian(std::vector<std::vector<Point>> dataSegemntedPoint, int window_size){
    std::vector<std::vector<Point>> dataSegemntedPoint_smooth;
    for(int i = 0; i < dataSegemntedPoint.size(); i++){
        std::vector<Point> segment_smooth;
        for(int j = 0; j < dataSegemntedPoint[i].size(); j++){
            if (j < window_size || j >= dataSegemntedPoint[i].size() - window_size) {
                // Se il dato è troppo vicino al bordo, non è possibile calcolare la media
                segment_smooth.push_back(dataSegemntedPoint[i][j]);
            }else {
                // Calcola la media dei dati nel range
                double total_y = 0;
                for (int k = j - window_size; k < j + window_size + 1; ++k) {
                    total_y += dataSegemntedPoint[i][k].y;
                }
                float mean = total_y / (window_size * 2 + 1);
                segment_smooth.push_back(Point(dataSegemntedPoint[i][j].x, mean, true));
            }
        }
        dataSegemntedPoint_smooth.push_back(segment_smooth);
    }
    return dataSegemntedPoint_smooth;
}

std::vector<std::vector<Point>>calc_first_derivative_segment_smooth(std::vector<std::vector<Point>> data){
    std::vector<std::vector<Point>> first_derivative_segment_smooth;
    for(int i = 0; i < data.size(); i++){
        std::vector<Point> derivative_segment_smooth;
        for(int j = 0; j < data[i].size()-1; j++){
            derivative_segment_smooth.push_back(Point(data[i][j].x,findiff(data[i][j+1].y, data[i][j].y, data[i][j+1].x - data[i][j].x), true));
        }
        first_derivative_segment_smooth.push_back(derivative_segment_smooth);
    }
    return first_derivative_segment_smooth;
}





std::vector<std::vector<Point>> circle_centers(std::vector<float> ranges, float incrementTheta, float thetaBegin){
    std::vector<float> ranges_clean = clean_first_20(ranges);


    //calcualte the first derivative of the polar data
    std::vector<float> first_derivative_polar;
    calc_first_derivative_polar(ranges_clean, first_derivative_polar,incrementTheta);

    //divide in segment
    std::vector<std::vector<float>> segmentedData;
    std::vector<int> segment_index_begin;

    divide_in_segment(ranges_clean,first_derivative_polar, segmentedData, segment_index_begin);

    //trasform each segment in cartesian
    std::vector<std::vector<Point>> segmentedDataCartesian;
    float const initial_theta_centerd = 3.14/2;// 90 degree
    segmentPolar_to_segmentCartesian(segmentedData, segmentedDataCartesian,initial_theta_centerd, incrementTheta);

    //smooth segenment
    std::vector<std::vector<Point>> segmentedDataCartesian_smooth;
    segmentedDataCartesian_smooth = smooth_data_cartesian(segmentedDataCartesian, WINDOW_SIZE_SMOOTH);


    //calcualte the derivative of each segment_smooth
    std::vector<std::vector<Point>> first_derivative_segment_smooth;
    first_derivative_segment_smooth = calc_first_derivative_segment_smooth(segmentedDataCartesian_smooth);

 
    //part to find the segment that are parabolas
    
    
    
    std::vector<int> segment_filtered_index;
    //chek is the fist derivative is neato to zero
    for(int i = 0; i < first_derivative_segment_smooth.size(); i++){
        int count = 0;
        for(int j = 0; j < first_derivative_segment_smooth[i].size(); j++){
            if(std::fabs(first_derivative_segment_smooth[i][j].y) < 0.1){
                count ++;
            }
        }
        if(count<10 && count > 0){ //the count threshold defiend from observation of the data plot
            segment_filtered_index.push_back(i);
        }
    }

    //if the segment are too short remove it
    std::vector<int> segment_filtered_index_long;
    for(int i = 0; i < segment_filtered_index.size(); i++){
        if(segmentedDataCartesian_smooth[segment_filtered_index[i]].size() > 6){ //TODO: the count threshold are hardocded
            segment_filtered_index_long.push_back(segment_filtered_index[i]);
        }
    }

    //if the derivative after the fist zero is positive add to the list
    std::vector<int> segment_filtered_index_positive;
    for(int i = 0; i < segment_filtered_index_long.size(); i++){
        float total = 0;
        int count = 0;
        for(int j = 0; j < first_derivative_segment_smooth[segment_filtered_index_long[i]].size(); j++){ //cheking form the beginnig beacuse the vactore is reverded
            if(!(std::fabs(first_derivative_segment_smooth[segment_filtered_index_long[i]][j].y) < 0.1)){
                total += first_derivative_segment_smooth[segment_filtered_index_long[i]][j].y;
                count ++;
            }
            else{
                break;
            }
        }
        double mean = total/count;
        if(mean > 0){
            segment_filtered_index_positive.push_back(segment_filtered_index_long[i]);
        }
    }

    //if the derivative from the end the fist zero(from the end) is negative add to the list
    std::vector<int> segment_filtered_index_negative;
    for(int i = 0; i < segment_filtered_index_positive.size(); i++){
        float total = 0;
        int count = 0;
        for(int j = first_derivative_segment_smooth[segment_filtered_index_positive[i]].size()-1; j > 0; j--){ //cheking form the beginnig beacuse the vactore is reverded
            if(!(std::fabs(first_derivative_segment_smooth[segment_filtered_index_positive[i]][j].y) < 0.1)){
                total += first_derivative_segment_smooth[segment_filtered_index_positive[i]][j].y;
                count ++;
            }
            else{
                break;
            }
        }
        double mean = total/count;
        if(mean < 0){
            segment_filtered_index_negative.push_back(segment_filtered_index_positive[i]);
        }
    }

    std::vector<int> index_circle_segment = segment_filtered_index_negative;

    //transofrm index referred to the segment index to the index referred to the original data
    std::vector<int> index_circle_segment_original,lengthSegment_circle_segment_original;
    for(int i = 0; i < index_circle_segment.size(); i++){
        index_circle_segment_original.push_back(segment_index_begin[index_circle_segment[i]]);
        lengthSegment_circle_segment_original.push_back(segmentedData[index_circle_segment[i]].size());
    }


    //crate polar vectors with circle segment
    std::vector<std::vector<float>> circle_segment_polar;
    for(int i = 0; i < index_circle_segment_original.size(); i++){
        std::vector<float> segment;
        for(int j = 0; j < lengthSegment_circle_segment_original[i]; j++){
            segment.push_back(ranges_clean[index_circle_segment_original[i]+j]);
        }
        circle_segment_polar.push_back(segment);
    }

    //polar to cartesian segment with the correct reference
    std::vector<std::vector<Point>> circle_segment_cartesian;
    for(int i = 0; i < circle_segment_polar.size(); i++){
        std::vector<Point> segment;
        float theta = thetaBegin + index_circle_segment_original[i]*incrementTheta;
        for(int j = 0; j < circle_segment_polar[i].size(); j++){
            segment.push_back(Point(circle_segment_polar[i][j], theta + (j*incrementTheta)));
        }
        circle_segment_cartesian.push_back(segment);
    }


    return circle_segment_cartesian;
}


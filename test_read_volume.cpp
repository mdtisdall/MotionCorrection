#include "Weighted_Gauss_Newton_Ref_Grad.h"
#include "Weighted_Gauss_Newton_New_Grad.h"

#include "Algorithms.h"

#include "BinaryFile.h"
#include "ReadVolume.h"

#include <iostream>
#include <fstream>

#include <string>
#include <sstream>


using namespace Algorithms;

int main(){
    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int cubeSize = 32;
    VolumeT vol(cubeSize); 
    int bytesRead = ReadVolume<VolumeT>::read_volume(&vol, path, cubeSize);
    
    // std::cout << "bytes Read = " << bytesRead << std::endl;
    /* print volume */
    // for (int i = 0; i < cubeSize*cubeSize*cubeSize; i++){
    //     std::cout << " " << vol.at(i); 
    // }

    string path1 = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/testVol1.dat";
    VolumeT vol_ref(cubeSize);
    int bytesRead_ref = BinaryFile<VolumeT>::read(&vol_ref, path1);
    for (int i = 0; i < vol_ref.size(); ++i){
        if (vol.at(i) != vol_ref.at(i)){
            cout << i << " ";
        }
    }
    std::cout << endl;
    return 0;
}

#ifndef ReadVolume_h
#define ReadVolume_h

#include "BinaryFile.h"


template <typename AtAddressableT>
class ReadVolume{
  public:
    typedef typename AtAddressableT::value_type value_type;
    typedef float dataT;
    typedef std::complex<dataT> complexT;
    typedef std::vector<complexT> vectorT;

    static int read_volume(AtAddressableT *buffer, std::string path, int cubeSize){
        vectorT readbuffer(cubeSize*cubeSize);
        int bytesRead = 0;
        for (int i = 0; i < cubeSize; ++i){
            std::stringstream ss;
            ss << i;
            std::string slice_path = path+"_slice_" + ss.str() +".dat";
            int inputFile = open(slice_path.c_str(), O_RDONLY);

            if(-1 == inputFile) {
                std::cerr << "Could not open file: " << slice_path << std::endl;
                exit(1);
            }

            bytesRead +=
                ::read(inputFile, &(readbuffer.at(0)), sizeof(complexT) * readbuffer.size());

            dataT *writebuffer = &(buffer->at(0)) + cubeSize*cubeSize*i;


            for(int j = 0; j < cubeSize*cubeSize; ++j){
                writebuffer[j] = sqrt(norm(readbuffer[j]));
            }
            close(inputFile);
        }
        return bytesRead/2;
    }

};
#endif

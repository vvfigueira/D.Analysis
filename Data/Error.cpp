#include <fstream>
#include "iostream"

void Change(char* name){

    std::fstream eFile;
    eFile.open(name, std::ios::in);
    std::string name1 = name;
    int n = name1.size();
    name1.erase(n-4,4);
    name1= name1+"F.tsv";
    char *name2 = &name1[0];
    std::fstream aFile (name2, std::fstream::app);
    double x[4];
    while(1){

        eFile >> x[0] >> x[1] >> x[2] >> x[3];
        aFile << x[0]/2 - x[2]/2 <<"\t" <<x[1]/2-x[3]/2 <<"\n";
        if(eFile.eof()) break; 
        
    }
    aFile.close();
    eFile.close();
}

int main(int argc, char** argv){

    switch (argc){
        case 2:
            Change(argv[1]);
            break;
        
        default:
            std::cout << "Modo de Uso: ./Error [Nome_ do_Arquivo]\n";
            break;
    }

    return 0;
}
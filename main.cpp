#include <iostream>
#include <cstdint>
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <vector>
#include <bsoncxx/json.hpp>
#include <mongocxx/client.hpp>
#include <mongocxx/uri.hpp>
#include <mongocxx/stdx.hpp>
#include <mongocxx/uri.hpp>
#include <mongocxx/instance.hpp>
#include <bsoncxx/builder/stream/helpers.hpp>
#include <bsoncxx/builder/stream/document.hpp>
#include <bsoncxx/builder/stream/array.hpp>
#include "parasail.h"
#include "parasail/matrices/blosum62.h"
#include "parasail/matrix_lookup.h"
#include <mongocxx/pool.hpp>
#include <stdio.h>
#include <omp.h>
//#include <mpi.h>
using bsoncxx::builder::stream::close_array;
using bsoncxx::builder::stream::close_document;
using bsoncxx::builder::stream::document;
using bsoncxx::builder::stream::finalize;
using bsoncxx::builder::stream::open_array;
using bsoncxx::builder::stream::open_document;


void parasailMtx(const char *inputProt1, const char *inputProt2, const char *nameOne, const char *nameTwo);



int main() {
    mongocxx::instance instance{};
    mongocxx::client client{mongocxx::uri{}};

    mongocxx::database db = client["AK"];
    mongocxx::collection coll = db["NCBI"];

    mongocxx::collection coll2 = db["EPI2"];

    mongocxx::cursor cursor = coll.find({});
    mongocxx::cursor cursor2 = coll2.find({});

    std::cout << typeid(cursor).name() << std::endl;

    std::vector<std::string> cursorList = {};
    std::vector<std::string> cursorNames = {};

    std::vector<std::string> cursorList2 = {};
    std::vector<std::string> cursorNames2 = {};

    /*Parallel computing. MPI MULTIPROCESSING*/

    /*int rank, comm_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);*/



    for(auto const& doc : cursor) {
        //std::cout << "TEST FIRST PROT" << bsoncxx::to_json(doc) << "\n";
        bsoncxx::document::element protOne = doc["FASTA"];
        std::size_t posOne = std::string (protOne.get_utf8().value).find("\n");
        std::string seqOne = std::string (protOne.get_utf8().value).substr(posOne);
        cursorList.push_back(seqOne);

        bsoncxx::document::element nameOne = doc["NCBI_ID"];
        std::string nameVal = std::string (nameOne.get_utf8().value);

        cursorNames.push_back(nameVal);

    }

    for(auto const& doc : cursor2){
        bsoncxx::document::element protTwo  = doc["Epitope_Antigen_Sequence"];
        std::string prot2 = std::string (protTwo.get_utf8().value);
        cursorList2.push_back(prot2);

        bsoncxx::document::element nameTwo = doc["_id"];
        std::string  nameVal2 = std::string(nameTwo.get_utf8().value);

        cursorNames2.push_back(nameVal2);
    }

    #pragma omp parallel for num_threads(40)
    for(int i = 0; i < cursorList.size(); i++)
    {
        for(int j = 0; j < cursorList2.size(); j++) {
            //if (i%comm_size != rank) continue;
            const char *valueOne = cursorList[i].c_str();
            const char *valueTwo = cursorList2[j].c_str();
            std::cout << "CKECKING PROTS: " << valueOne << " AND " << valueTwo << "\n";
            const char *nameValOne = cursorNames[i].c_str();
            const char *nameValTwo = cursorNames2[j].c_str();

            std::string fNamePt1 = nameValOne;
            std::string fNamePt2 = nameValTwo;

            std::string fileName = "/home" + fNamePt1 + "_" + fNamePt2 + ".txt";
            const char *myFile = fileName.c_str();

            std::ifstream myFil(myFile);
            if (myFil.fail()) {
                //File does not exist code here
                std::cout << "THE COMPARED ARE: " << nameValOne << "    AND    " << nameValTwo << "\n";
                parasailMtx(valueOne, valueTwo, nameValOne, nameValTwo);
            }
        }
    }

    /*    for (auto const& doc2 : cursor) {
            std::cout << "FIRST PROTEIN" << bsoncxx::to_json(doc) << "\n";
            std::cout << "SECOND PROTEIN" << bsoncxx::to_json(doc2) << "\n";
            bsoncxx::document::element protOne = doc["FASTA"];
            //std::cout << "PROT ONE" << protOne.get_utf8().value << "\n";
            bsoncxx::document::element protTwo = doc2["FASTA"];
            //std::cout << "PROT TWO" << protTwo.get_utf8().value << "\n";

            //Test conv to const char:

            std::size_t posOne = std::string (protOne.get_utf8().value).find("\n");
            std::size_t posTwo = std::string (protTwo.get_utf8().value).find("\n");
            std::string seqOne = std::string (protOne.get_utf8().value).substr(posOne);
            std::string seqTwo = std::string (protTwo.get_utf8().value).substr(posTwo);
            //std::cout<<"SEQ DATA"<<"\n";
            //std::cout<<seqOne<<"\n";
            //std::cout<<seqTwo<<"\n";
            const char* valueOne = seqOne.c_str();
            const char* valueTwo = seqTwo.c_str();
            //std::cout<<"CKECKING PROTS: "<<valueOne<<" AND "<<valueTwo<<"\n";
            //Get names (IDs) of the compared proteins

            bsoncxx::document::element nameOne = doc["NCBI_ID"];
            bsoncxx::document::element nameTwo = doc2["NCBI_ID"];

            //Convert to text

            const char * nameValOne = std::string(nameOne.get_utf8().value).c_str();
            const char * nameValTwo = std::string(nameTwo.get_utf8().value).c_str();


            parasailMtx(valueOne, valueTwo, nameValOne, nameValTwo);

            //parasailMtx();

        }
    }


    */


    //std::cout << "Finished!" << std::endl;
    //MPI_Finalize();
    return 0;
}
void parasailMtx(const char *inputProt1, const char *inputProt2, const char *nameOne, const char *nameTwo){
    //std::cout<<"FUN WRKS"<<"\n";
    //std::cout<<"ACTUAL LEN"<< strlen(inputProt1) << "\n";
    std::cout<<"ACTUAL PROT"<<inputProt1 << "     AND ANOTHER        " << inputProt2<<"\n";
    const char *s1 = inputProt1;
    const char *s2 = inputProt2;
    int s1Len = (int)strlen(s1);
    int s2Len = (int)strlen(s2);
    parasail_result_t *result;
    const parasail_matrix_t *matrix = NULL;


    /* note the address-of operator '&' */
    //result = parasail_sw(s1, s1Len, s2, s2Len, 11, 1, &parasail_blosum62);
    //parasail_result_free(result);
    matrix = parasail_matrix_lookup("blosum62");
    result = parasail_sg_trace_striped_sat(
            s1, s1Len, s2, s2Len, 11, 0, matrix); //)//parasail_nw_trace_striped_sse41_128_64(

    //std::cout << typeid(result).name() << std::endl;
    //std::cout<<result<<"\n";


    //result = parasail_sw(s1, s1Len, s2, s2Len, 11, 1, matrix);
    //parasail_result_free(result);
    //Tests length
    //std::cout<<"LENGTH"<<"\n";
    //std::cout<<s1Len<<"\n";
    //std::cout<<s2Len<<"\n";
    std::string fNamePt1 = nameOne;
    std::string fNamePt2 = nameTwo;
    std::string fileName = "/home/" + fNamePt1 + "_" + fNamePt2 + ".txt";
    const char *myFile = fileName.c_str();
    std::cout<< " THE NEW FILE NAME "<< fileName;
    FILE *myfile = fopen(myFile, "w+");

    if (s1Len > s2Len) {
        std::cout<<"TEST IF"<<"\n";
        std::cout<<"Preliminary result is"<<result<<"\n";
        std::cout<<"CHECKING THE LEN VAL   "<<s1Len<<"\n";
        std::cout<<"Input len 1 \n"<<s1Len<<"\n ANd 2 len: \n"<<s2Len<<"\n";
        parasail_traceback_generic_extra(
                s1, s1Len, s2, s2Len,
                nameOne, nameTwo, matrix,
                result, '|', '+', '-', s1Len, 20, 1, s1Len, myfile);
        std::cout<<"TEST AFTER"<<"\n";
        //fclose(myfile);

    }
    else {
        std::cout<<"TEST ELSE"<<"\n";
        parasail_traceback_generic_extra(
                s1, s1Len, s2, s2Len,
                nameOne, nameTwo, matrix,
                result, '|', '+', '-', s2Len, 20, 1, s2Len, myfile);
        std::cout<<"TEST AFTER ELSE"<<"\n";
        //fclose(myfile);

    }
    parasail_result_free(result);
    //parasail_matrix_free(*matrix);
    //std::cout<<"LEFT ELSE"<<"\n";
    //std::cout<<result<<"\n";
    fclose(myfile);
}

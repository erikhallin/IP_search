#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <windows.h>

using namespace std;

struct st_ip
{
    st_ip()
    {
        res='X';
        pos=-1;
        ip=-1;
    }
    st_ip(char _res,int _pos,float _ip)
    {
        res=_res;
        pos=_pos;
        ip=_ip;
    }

    char res;
    int  pos;
    float ip;
};

float get_ip(string seq);

int main()
{
    cout<<"pI Search Software\n\nCalculate the theoretical pI of a window of a protein sequence.\n\n";

    string protein_seq;
    string line;
    int run_mode=0;//1=Ntrunc,2=Ctrunc,3=window
    int window_size=0;

    //read input seq
    ifstream input_file("input_seq.txt");
    if(input_file==0)
    {
        cout<<"ERROR: A protein sequence must be present in a file named input_seq.txt\n";
        return 1;
    }
    while(getline(input_file,line))
    {
        if(line[0]=='>') continue;
        protein_seq.append(line);
    }
    input_file.close();
    //go through seq
    for(int i=0;i<(int)protein_seq.length();i++)
    {
        protein_seq[i]=toupper(protein_seq[i]);
    }

    cout<<"Protein sequence input:\n"<<protein_seq<<endl;

    //ask search options
    cout<<"\nSelect run mode, (1) N-truncations, (2) C-truncations, (3) Window search: ";
    string user_input;
    getline(cin,line);
    switch(line[0])
    {
        case '1': run_mode=1; break;
        case '2': run_mode=2; break;
        case '3': run_mode=3; break;
    }
    if(run_mode==3)
    {
        cout<<"Select search window size: ";
        getline(cin,line);
        window_size=(int)atof(line.c_str());

        if(window_size<=0)
        {
            cout<<"ERROR: Bad window size\n";
            return 1;
        }
    }
    if(run_mode!=1&&run_mode!=2&&run_mode!=3)
    {
        cout<<"ERROR: Bad selection\n";
        return 1;
    }

    //scan sequence
    vector<st_ip> vec_ips;
    for(int i=0;i<(int)protein_seq.length();i++)
    {
        switch(run_mode)
        {
            case 1://ntrunc
            {
                string sub_seq(protein_seq,i);
                vec_ips.push_back(st_ip(protein_seq[i],i+1,get_ip(sub_seq)));

                //cout<<sub_seq<<endl;
            }break;

            case 2://ctrunc
            {
                //string sub_seq(protein_seq,0,(int)protein_seq.length()-i);
                //vec_ips.push_back(st_ip(protein_seq[i],i+1,get_ip(sub_seq)));

                string sub_seq(protein_seq,0,i);
                vec_ips.push_back(st_ip(protein_seq[i],i+1,get_ip(sub_seq)));

                //cout<<sub_seq<<endl;
            }break;

            case 3://window
            {
                string sub_seq;
                int start_pos=i-window_size;
                if(start_pos<0) start_pos=0;
                for(int j=start_pos;j<(int)protein_seq.length()&&j<=i+window_size;j++)
                {
                    sub_seq.append(1,protein_seq[j]);
                }
                vec_ips.push_back(st_ip(protein_seq[i],i+1,get_ip(sub_seq)));

                //cout<<sub_seq<<endl;
            }break;
        }
    }

    //print results
    ofstream output_file("output.txt");
    if(output_file==0)
    {
        cout<<"ERROR: Could not create output file\n";
        return 1;
    }
    output_file<<"Number\tResidue\tpI"<<endl;
    for(int i=0;i<(int)vec_ips.size();i++)
    {
        output_file<<vec_ips[i].pos<<"\t"<<vec_ips[i].res<<"\t"<<vec_ips[i].ip<<endl;
    }
    output_file.close();

    cout<<"\nComplete, results printed in output.txt\n\n";

    cout<<"Plot results? (Y/N): ";
    getline(cin,line);
    if(line[0]=='y'||line[0]=='Y')
    {
        ofstream plot_file("plot.data");
        if(plot_file==0)
        {
            cout<<"ERROR: Could not create a file\n";
            return 1;
        }
        plot_file<<"Number\tpI"<<endl;
        for(int i=0;i<(int)vec_ips.size();i++)
        {
            plot_file<<vec_ips[i].pos<<"\t"<<vec_ips[i].ip<<"\t"<<0<<endl;
        }
        plot_file.close();

        system("graph_window plot.data");

        remove("plot.data");
    }

    return 0;
}

float get_ip(string seq)
{
    char Asp = 'D';
    char Glu = 'E';
    char Cys = 'C';
    char Tyr = 'Y';
    char His = 'H';
    char Lys = 'K';
    char Arg = 'R';

    int AspNumber = 0;
    int GluNumber = 0;
    int CysNumber = 0;
    int TyrNumber = 0;
    int HisNumber = 0;
    int LysNumber = 0;
    int ArgNumber = 0;

    for (int i = 0; i <=(int)seq.length() - 1; i++)
    {
        if (seq[i] == Asp) AspNumber++;
        if (seq[i] == Glu) GluNumber++;
        if (seq[i] == Cys) CysNumber++;
        if (seq[i] == Tyr) TyrNumber++;
        if (seq[i] == His) HisNumber++;
        if (seq[i] == Lys) LysNumber++;
        if (seq[i] == Arg) ArgNumber++;
    }

    float NQ = 0.0; //net charge in given pH
    float QN1=0;  //C-terminal charge
    float QN2=0;  //D charge
    float QN3=0;  //E charge
    float QN4=0;  //C charge
    float QN5=0;  //Y charge
    float QP1=0;  //H charge
    float QP2=0;  //NH2 charge
    float QP3=0;  //K charge
    float QP4=0;  //R charge

    float pH = 0.0;

    while(true)
    {
        QN1=-1/(1+pow(10,(3.65-pH)));
        QN2=-AspNumber/(1+pow(10,(3.9-pH)));
        QN3=-GluNumber/(1+pow(10,(4.07-pH)));
        QN4=-CysNumber/(1+pow(10,(8.18-pH)));
        QN5=-TyrNumber/(1+pow(10,(10.46-pH)));
        QP1=HisNumber/(1+pow(10,(pH-6.04)));
        QP2=1/(1+pow(10,(pH-8.2)));
        QP3=LysNumber/(1+pow(10,(pH-10.54)));
        QP4=ArgNumber/(1+pow(10,(pH-12.48)));

        NQ=QN1+QN2+QN3+QN4+QN5+QP1+QP2+QP3+QP4;

        if (pH>=14.0)
        {                                                 //you should never see this message
            cout<<"Warning: pH is higher than 14"<<endl;  //
            break;
        }
        if (NQ<=0)                            // if this condition is true we can stop calculate
        break;

        pH+=0.01;                            // if not increase pH
    }

    return pH;
}

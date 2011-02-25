//------------------DataDumper--------------------//

#ifndef _DATA_DUMPER_HPP_
#define _DATA_DUMPER_HPP_

#include <omd/detector.hpp>

using std::string;

#define EXLN 5

/**
 * @ingroup detector
 * @brief The detector to write data to file
 * 
 * Detectors class which prints out the position and value of Atoms
 * to a file. (not to be confused with SysMonitor) One can define 
 * the main name of data files in the constructor.
 * The extention name of these files is the dumping sequence. One can also 
 * define the format the width of extention such as "000" or "0000" by giving
 * the default parameter "extlen" the number of digit he want. Giving 0 (zero) 
 * to this parameter will cause this class to use only bare numbers in the file
 * extentions (e.g. 1, 2, 3, .... instead of 001, 002, 003, ... etc).
 *
 * Parameters:
 * - dump.sample (t) : sample time (in ps).
 * - dump.filename (file name) : the data filename without extension.
 * - dump.xlength (n) : the length of enumerated name extension.
 *
 * The default sample time is the initial time step. The default will also be used
 * if the sample time is less than the initial time step.
 *
 * Caution! By default this data dumper will dump the first data file in
 * the file name with 001 extension.
 *
 * To explain:
 *  - main detector
 *  - joining detectors, why?
 *  - file naming
 *
 *
*/

class DataDumper: public Detector {
protected:
    OMD_INT FileSeq, FileMax;
    OMD_INT XLength;
    string FieldsDef;
    DataDumper* MainDumper;

public:
    
   	DataDumper(string TargetAtom="",OMD_FLOAT tm=-1.0,string fn="Data",OMD_INT extlen=3):
   	Detector(tm) {
        TargetName=TargetAtom;
		SetFilename(fn);
        FileSeq=0;
        XLength=(extlen<EXLN)?extlen:EXLN;
        if (XLength==0) FileMax=100000;
        else FileMax=OMD_INT(pow(10, XLength+1));
        MainDumper=NULL;
        
        set_name("DATA DUMPER");
        register_class(get_name());
        
    }
    
    /** 
     * Join is used to define the main data dumper.
     * By joining to other data dumper, the class will not take responsible
     * to write the output file. Only the ones which have MainDumper==NULL, will
     * create output files. In order to use join capability, the data must be 
     * store in the aux variable of atoms.
     * 
     */

    Detector* Join(MDGadget* det){

    	assert(det, "Can not find DataDumper to join");
    	assert(det->type_of("data_dumper"), "Gadget "+det->get_name()+
    	       " is not a data dumper. can not join.");
    	       
    	MainDumper=dynamic_cast<DataDumper*>(det);
    	return this;
	}
 
    virtual string GetFilename() {
        OMD_INT i=0,j;
        OMD_CHAR st[64];
		
		if (XLength==0)	return (FilenamePrefix+Filename);
			
        SetFilenamePostfix("00000");
        FileSeq+=1;
        sprintf(st, "%d", FileSeq);
        if (FileSeq<FileMax) {
	        if (XLength==0) SetFilenamePostfix(st);
		    else {
            	while (st[i]!=0x0) i++;
            	for (j=1;j<=i;j++) FilenamePostfix[XLength-j]=st[i-j];
            	FilenamePostfix[XLength]=0x0;
        	}
        } else SetFilenamePostfix("---");
        
    	return (FilenamePrefix+Filename+"."+FilenamePostfix);
	}
	
	virtual void Init(MDSystem* WorkSys) {
		Detector::Init(WorkSys);

		SysParam->peek("dump.sample", TSample);
		SysParam->peek("dump.filename", Filename);
		SysParam->peek("dump.xlength", XLength);
		
		if(MainDumper!=NULL){
			assert(Target==MainDumper->Target,
			       "joining data dumpers with with different target");
			TSample=MainDumper->GetSampleTime();
			NextSample=TSample;
			Target->SetWriteMode(WM_VELOCITY);
		}
		
		if(MainDumper==NULL) RestartVariable("file_sequence", FileSeq);

	}

	virtual void PrintInfo(ostream& ost) {
		Detector::PrintInfo(ost);
		ost << "  sample_time="<<TSample<<"; filename="<<Filename<<"; extension_length="<<XLength<<std::endl;
	}

    virtual bool Check(){
    	if(MainDumper==NULL){
			for(OMD_INT a=0;a<MAXAUXVAR;a++)
				if(System->PrintableAux[a])
					FieldsDef.append(" "+System->AuxNameTag.at(a));
		}
		return true;
	}

    virtual void Measure() {
    	if(MainDumper==NULL){
    		Target->PushInfo("$ Time "+as_string(System->ElapsedTime,"%0.5e"));
        	Target->DumpAtoms(GetFilename(),0,
        		System->PrintableAux,
        		System->AuxFormat,
        		FieldsDef);
        	Target->PopInfo();
		}
    }
};

#endif

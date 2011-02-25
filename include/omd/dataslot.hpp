/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 3.0 (2009,2011) 
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 * ObjectMD header file
 *
 *      Screen data printing manager
 *
*/

#ifndef _DATA_SLOT_H_
#define _DATA_SLOT_H_

/**
 * @ingroup baseclass
 * @brief Data Slot
 *
 * The class is used to register a data entity to the system. MDSystem has two
 * data slot arrays, MessageSlots and RestartVars. The first is used to display
 * data to the screen, and the second registers the data that must be saved in
 * the restart binary file.
 * 
 * A gadget must register data to be displayed using
 * MDGadget::RegisterMsgSlot() function. If a data needed to be saved in the
 * binary file, MDGadget::RestartVariable() function must be invoked. This 
 * function will set the variable pointer the the one given in the parameter, 
 * and if a slot is already exist, the data will be assigned as well.
 * 
 */

class DataSlot {

	string *StrData;
	OMD_FLOAT *DblData;
	OMD_INT    *IntData;
	OMD_INT priority;

	string default_data;
	string Label;
	string Format;
	enum {m_undefined, m_int, m_float, m_string} DataType;
	
	OMD_CHAR   FormattedText[512];
	bool   Printable;

public:
	
	/**
	 * label=text label to display. priority=order, lower means higher priority.
	 */
	
	DataSlot(string vlabel="", OMD_INT vpriority=1) {
		Label.assign(vlabel);
		priority=vpriority;
		FormattedText[0]='\0';
		Printable=true;
		DataType=m_undefined;
	}
	
	virtual ~DataSlot(){}	
	
	OMD_INT GetPriority(){return priority;}
	
	bool    IsPrintable(){return Printable;}
	
	OMD_INT AsInt(){
		switch(DataType) {
			case m_int:
				return *IntData;
			case m_float:
				return (OMD_INT)(*DblData);
			case m_string:
				return(atoi(StrData->c_str()));
			case m_undefined:
				throw "undefined message type!";
		}
		return 0; // avoids warning
	}
		
	OMD_FLOAT  AsDouble(){
		switch(DataType) {
			case m_int:
				return (OMD_FLOAT)(*IntData);
			case m_float:
				return *DblData;
			case m_string:
				return(atof(StrData->c_str()));
			case m_undefined:
				throw "undefined message type!";
		}
		return 0.0; // avoids warning
	}
	
	string AsString(){
		OMD_CHAR st[128];
		string retst;

		switch(DataType) {
			case m_int:
				if(Format.empty())sprintf(st, "%d", *IntData);
				else sprintf(st, Format.c_str(), *IntData);
				retst.assign(st);
				break;
			case m_float:
				if(Format.empty())sprintf(st, "%f", *DblData);
				else sprintf(st, Format.c_str(),*DblData);
				retst.assign(st);
				break;
			case m_string:
				retst.assign(*StrData);
				break;
			case m_undefined:
				throw "undefined message type!";
			}

		return retst;
	}
	
	void* GetDataPointer() {
		switch(DataType) {
			case m_int:
				return IntData;
			case m_float:
				return DblData;
			case m_string:
				return StrData;
			case m_undefined:
				throw "undefined message type!";
		}
	}
	
	DataSlot* SetFormat(const OMD_CHAR* fmt){Format.assign(fmt);return this;}
	DataSlot* SetPrintable(bool p){Printable=p;return this;}
	DataSlot* SetData(OMD_SIZET& dat){IntData=(OMD_INT*)&dat;DataType=m_int;return this;}
	DataSlot* SetData(OMD_INT& dat) {IntData=&dat;DataType=m_int;return this;}
	DataSlot* SetData(OMD_FLOAT& dat) {DblData=&dat;DataType=m_float;return this;}
	DataSlot* SetData(string& dat){StrData=&dat;DataType=m_string;return this;}
	DataSlot* SetLabel(const OMD_CHAR* dat){Label.assign(dat);return this;}
	
	DataSlot* SetDefaultData(string dat){
		default_data=dat;
		SetData(default_data);
		return this;
	}

	string& GetLabel(){return Label;}
	OMD_CHAR* GetFormattedText() {
		string st;
		if(Label.empty())st.assign(AsString());
		else st.assign(Label+":"+AsString());		
		memset(FormattedText,0,512);
		st.copy(FormattedText,511);
		return FormattedText;
	}
};

#endif

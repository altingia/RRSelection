#ifndef calculate_H_
#define calculate_H_


using namespace std;


/////////////////////////                 method1 /////////////////////////////////
int cal_RR_MA ( vector<BaseType>  & Base1  , vector<BaseType> &  Base2 ,  double &  CalResult, statementVar & Var )
{
	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;
	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++ ;
	}
	
	Var.tmpAA=(Var.DDE[1][1])+(Var.DDE[1][0]);
	if (Var.tmpAA==0)
	{
		return 0 ;
	}
	if ( (Var.DDE[1][1]+Var.DDE[0][1])==0)
	{
		return 0 ;
	}
	
	Var.ALL_count=Var.DDE[0][0]+Var.DDE[0][1]+Var.tmpAA;
	Var.probHaps[0]=((Var.DDE[0][0])/Var.ALL_count);
	Var.probHaps[1]=((Var.DDE[0][1])/Var.ALL_count);
	Var.probHaps[2]=((Var.DDE[1][0])/Var.ALL_count);

	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];

	Var.Cal_B=(Var.pA1)*(Var.pA2);
	Var.Cal_A  = 1.0-(Var.pA1+Var.pA2)+Var.Cal_B ;

	if  (Var.Cal_A==0  || Var.Cal_B==0 )
	{
		if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
		if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
		if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}

		Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
		Var.pA2 = Var.probHaps[0]+Var.probHaps[2];

		Var.Cal_B = (Var.pA1)*(Var.pA2);
	    Var.Cal_A  = 1.0-(Var.pA1+Var.pA2)+Var.Cal_B ;

	}

	Var.D_A = Var.probHaps[0]-Var.Cal_B ;
	CalResult = (Var.D_A*Var.D_A)/(Var.Cal_A*Var.Cal_B);

	return 1;
}





void PairWiseRRCal(In3str1v * paraFA04, Para_18 * para_18 , map <llong, vector <BaseType>  >  & SNPList ,  StarRsult  *  All_Stat[])
{
	double  CalResult;
	llong Dis=0;
	statementVar  Var;
	map<llong, vector <BaseType> >  :: iterator key2;
	map<llong, vector <BaseType> >  :: iterator key2_se;

	key2=SNPList.begin();
	Var.Asize= (key2->second).size();
	for ( ; key2!=SNPList.end(); key2++)
	{
		key2_se=key2 ; key2_se++ ;
		for( ; key2_se!=SNPList.end(); key2_se++)
		{
			Dis=(key2_se->first)-(key2->first);
			if ( Dis> (paraFA04->WindowSize))
			{
				break ;
			}
			if (cal_RR_MA( key2->second , key2_se->second  ,CalResult, Var) ==1)
			{
				int count=int((key2->first)/(paraFA04->bin));
				int SeCount=int(Dis/(paraFA04->bin));
				All_Stat[SeCount][count].Count++;
				All_Stat[SeCount][count].sumRR+=CalResult;
			}
		}
	}
}





#endif // calculate_H_  ;







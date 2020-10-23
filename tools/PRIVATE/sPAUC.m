function result = sPAUC(pAUC,FPRL,FPRU)
%SClaculates McCIsh standardized partial AUC
Amx=FPRU-FPRL;
Amn=(FPRU-FPRL).*(FPRU+FPRL)/.2;
result=(1+(pAUC-Amn)./(Amx-Amn))./2;
end

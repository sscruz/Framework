#include <skimmer.C>

void runSkim(TString name="", TString path="")
{
  if (path=="") return;
  skimmer* a = new skimmer(name,path);
  delete a;

}

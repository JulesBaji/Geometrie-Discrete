#include <DGtal/base/Common.h>

int main(int argc, char** argv)
{

  DGtal::trace.info() << "Hello from DGtal ";
  DGtal::trace.emphase() << "(version "<< DGTAL_VERSION << ")"<< std::endl;

  return 0;
}
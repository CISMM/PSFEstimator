/**
 * This file is where the program starts. The
 * main() function does very little other than
 * instantiate a QApplication and the application
 * main window.
 */

#include <qapplication.h>

#include "PSFEstimator.h"

int main(int argc, char** argv) {

  // Qt initialization.
  QApplication app(argc, argv);
  
  PSFEstimator mainWindow;
  mainWindow.show();
  
  return app.exec();
}

#include <QApplication>
#include "main_window.h"
#include <QTimer>

int main(int argc, char* argv[])
{
    QApplication application(argc, argv);
    QApplication::setApplicationName("Navier");
    MainWindow window;
    window.show();

    // Exercise widget construction and the event loop without user interaction.
    if (application.arguments().contains("--smoke-test")) {
        QTimer::singleShot(0, &application, &QApplication::quit);
    }
    return application.exec();
}

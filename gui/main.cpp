#include <QApplication>
#include <QLabel>
#include <QMainWindow>
#include <QTimer>

int main(int argc, char* argv[])
{
    QApplication application(argc, argv);
    QApplication::setApplicationName("Navier");
    QMainWindow window;
    window.setWindowTitle("Navier — Transient Model");
    auto* welcome = new QLabel("Navier\nTransient modeling");
    welcome->setAlignment(Qt::AlignCenter);
    window.setCentralWidget(welcome);
    window.resize(720, 480);
    window.show();

    // Exercise widget construction and the event loop without user interaction.
    if (application.arguments().contains("--smoke-test")) {
        QTimer::singleShot(0, &application, &QApplication::quit);
    }
    return application.exec();
}

#include <iostream>

#include "gui/Viewer.h"

#include <nsessentials/util/GLDebug.h>

int main()
{
	try
	{
#if _WIN32
		std::cout.imbue(std::locale("en-US"));
#else	
		std::cout.imbue(std::locale("en_US.utf8"));
#endif
	}
	catch (...)
	{
		std::cerr << "Warning: Could not set english locale." << std::endl;
	}

	nanogui::init();

	{
		nanogui::ref<Viewer> viewer = new Viewer();
		viewer->setVisible(true);

		nse::util::GLDebug::SetupDebugCallback();
		nse::util::GLDebug::IgnoreGLError(131185);

		try
		{
			nanogui::mainloop();
		}
		catch (std::runtime_error& e)
		{
			std::cerr << e.what() << std::endl;
		}

	}

	nanogui::shutdown();

	return 0;
}
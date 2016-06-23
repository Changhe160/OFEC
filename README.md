# OFEC
OPEN FRAMEWORK OF EVOLUTIONARY COMPUTATION

This is an open source project with the aim of providing an open framework for evolutionary computation (OFEC). OFEC includes fundamental structures of population based algorithms, problems, and algorithms' performance metrics. Based on this framework, you can easily implement your own algorithms and compare their performance with other algorithms. This project is totally for research purpose. I would like to encourage you to take part in this project and let's together make the research life easy, comfortable, and efficient for you and for everyone.

OFEC was implemented in C++11 language based on some libraries of Boost and several modern C++ design techniques. It supports solution representation in continuous and combinatorial space, or using hybird encoding schemes. It can be run in both Windows and Linux environments with a single thread or multi-thread mode.

OFEC has three basic components: algorithm, problem, and performance measurements. Several algorithm examples are provided, including three multi-objective optimization algorithms, the standard particle swarm optimization, the basic differential evolution, the original ant colony optimization, the standard genetic algorithm, and several algorithms of mine. The problem component collects problems mainly from continuous domain, including two dynamic optimization problems (the moving peak benchmark and the GDBG benchmark), global optimization problems, multi-modal optimization problems, and travelling salesman problems. Traditional performance metrics for dynamic optimization and global optimization are implemented in the performance component.

The current version of OFEC is 0.4.2 develped with MS Visual Studio Community 2015 and the BOOST version is 1.59, a makefile is also available for Linux users.

Please feel free to contact me by changhe.lw@gmail.com if you have any questions.

# Week notes

## Week 1

### Lecture 1, week 1, Monday

- Pick a set of class meeting times

- Sign up for the class piazza page using the piazza sidebar link on [our canvas site](https://canvas.ubc.ca/courses/180886)

- Review each of the links in the [course website background information](https://phaustin.github.io/a405/)

- try an install of the conda environment on your laptop following the instructions at [this page](https://phaustin.github.io/a405/python_install.html#laptop-installs)

- Read [Stull chapter 1 pp. 6-16](https://www.eoas.ubc.ca/books/Practical_Meteorology/prmet102/Ch01-atmos-v102b.pdf) and skim  Lohmann chapter 1. Think about how you would approach the Lohmann Chapter 1 long problems 4, 5, and 6 at the end of the chapter.  We will do these three problems in class together this week.


#### Learning objectives for Lohmann and Stull chapters 1

  - Identify principal cloud types
  - Explain how clouds can both heat and cool the planet
  - Build your intuition about droplet size and spacing in a typical cloud
  - Introduce hydrostatic balance and use it to determine the mass of the atmosphere
  
### Lecture 2, week 1, Thursday

  - Finish the python install, using conda-lock if conda is having trouble solving for the environment.yml file
  - Introduce my take on {ref}`sec:hydrostat`
  - Do the {ref}`worksheet1`

####  Do for Tuesday
  
Read the Thompkins text through p. 14 and start {ref}`assignment:adiabatic`

#### Learning objectives for the Thompkins reading

  - Distinguish between processes (heat and work) and state variables (energy, temperature, pressure, entropy etc.).  One way to think about this:  in thermodynamics, heat and work are always verbs ("I heat the burrito", "I work the problem") and never nouns ("the heat is killing me", "that's a lot of work").  State variables are never verbs ("don't pressure me"), they are always nouns ("the pressure is 850 hPa")
  - Conserved variables like mass, energy, and (for adiabatic processes) entropy are very powerful in physics, because they allow you to apply a constraint and simplify a problem by adding information.  Because of this, there is a whole zoo of thermodynamic variables that are conserved under various processes (precipitating, non-precipitating, freezing, etc.)
  
### Lecture 3, week 2, Tuesday

- Announcement:  I changed the Assignment 1 due date to Friday 11:30pm so there was time to ask questions during Thursday's class.

- Equation sheet:  For reference here's a pdf of the equation sheet that was staped to the back of the final exam in 2024:  [2024 equation sheet](https://www.dropbox.com/scl/fi/5scb458wnoymr4s8ta4vt/equation_sheet.pdf?rlkey=29ymawiljakqnv3wquqa4wm8e&st=4mmzb2ii&dl=0)

Learning objectives:

- Review the  Thompson reading through page 14.  I've added links to three more in-depth notes on this material:

    - [Kinetic temperature note](https://drive.google.com/file/d/1dBxcwxe0quNjYAihiFwMpF8kVmjmYh1p/view?usp=sharing).  This uses material from Physics 203 to define the temperature of an ideal gas in terms of the kinetic energy of the gas.
    - [First law notes](https://drive.google.com/file/d/12xontZRhRTKBjtG9y27RS6CH35eARcgp/view?usp=sharing).  This is a more detailed discussion of the how you can integrate the first law to get the internal energy and the entropy.
    - [Entropy and potential temperature](https://drive.google.com/file/d/1CjXm0AHYMJYmxf-Mm6CuWD7tIM_KBRI0/view?usp=sharing). This is a review of Thompson p. 14

2. Today's worksheet is about the A405 thermodynamic library:  [worksheet2.ipynb](https://drive.google.com/file/d/1NzrJy3mG9O9CGORxv3Vtujvd9CAQ-jwZ/view?usp=sharing).  There are 2 questions: one for Tuesday and one for Thursday.  You can hand in just once on Thursday.

For  Thursday:  Finish Assignment 1 and read my three notes.  Think about how you would actually get numbers for the internal energy, enthalpy and entropy of air in practice.  What are the limits of integration for temperature?  Are you starting at absolute zero?   How could that work if air would be rock solid at absolute zero?  Are the heat capacities really constant over that enormous temperature range?



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
  - Do the {ref}`worksheet1_solution`

####  Do for Tuesday
  
Read the Thompkins text through p. 14 and start {ref}`adiabatic_box_solution`

#### Learning objectives for the Thompkins reading

  - Distinguish between processes (heat and work) and state variables (energy, temperature, pressure, entropy etc.).  One way to think about this:  in thermodynamics, heat and work are always verbs ("I heat the burrito", "I work the problem") and never nouns ("the heat is killing me", "that's a lot of work").  State variables are never verbs ("don't pressure me"), they are always nouns ("the pressure is 850 hPa")
  - Conserved variables like mass, energy, and (for adiabatic processes) entropy are very powerful in physics, because they allow you to apply a constraint and simplify a problem by adding information.  Because of this, there is a whole zoo of thermodynamic variables that are conserved under various processes (precipitating, non-precipitating, freezing, etc.)
  
## Week 2

### Lecture 3, week 2, Tuesday

- Announcement:  I changed the Assignment 1 due date to Friday 11:30pm so there was time to ask questions during Thursday's class.

- Equation sheet:  For reference here's a pdf of the equation sheet that was staped to the back of the final exam in 2024:  [2024 equation sheet](https://www.dropbox.com/scl/fi/5scb458wnoymr4s8ta4vt/equation_sheet.pdf?rlkey=29ymawiljakqnv3wquqa4wm8e&st=4mmzb2ii&dl=0)

Learning objectives:

- Review the  Thompson reading through page 14.  I've added links to three more in-depth notes on this material:

    - [Kinetic temperature note](https://drive.google.com/file/d/1dBxcwxe0quNjYAihiFwMpF8kVmjmYh1p/view?usp=sharing).  This uses material from Physics 203 to define the temperature of an ideal gas in terms of the kinetic energy of the gas.
    - [First law notes](https://drive.google.com/file/d/12xontZRhRTKBjtG9y27RS6CH35eARcgp/view?usp=sharing).  This is a more detailed discussion of the how you can integrate the first law to get the internal energy and the entropy.
    - [Entropy and potential temperature](https://drive.google.com/file/d/1CjXm0AHYMJYmxf-Mm6CuWD7tIM_KBRI0/view?usp=sharing). This is a review of Thompson p. 14

2. Today's worksheet is about the A405 thermodynamic library:  [worksheet2.ipynb](https://drive.google.com/file/d/1NzrJy3mG9O9CGORxv3Vtujvd9CAQ-jwZ/view?usp=sharing).  There are 2 questions: one for Tuesday and one for Thursday.  You can hand in just once on Thursday.

#### For  Thursday:

Finish Assignment 1 and read my three notes.  Think about how you would actually get numbers for the internal energy, enthalpy and entropy of air in practice.  What are the limits of integration for temperature?  Are you starting at absolute zero?   How could that work if air would be rock solid at absolute zero?  Are the heat capacities really constant over that enormous temperature range?


### Lecture 4, week 2, Thursday

Learning objectives

1) Finish the second problem in {ref}`worksheet2_solution`

2) Start themodynamic diagram section 1.10 in the Thompson text using the {ref}`skew_coords` notebook.

3) Review derivatives and integrals in python with {ref}`derivs_ints`

4) Introduce {ref}`assign2_hydro_solution`

## Week 3

### Lecture 5, week 3, Tuesday

Learning objectives

1) Do the entropy calculation in the {ref}`worksheet2_solution` and review my notes on 
[entropy and potential temperature](https://drive.google.com/file/d/1CjXm0AHYMJYmxf-Mm6CuWD7tIM_KBRI0/view?usp=sharing)

2) I've turned the {ref}`derivs_ints` into today's worksheet.

3) Go over the {ref}`adiabatic_box_solution`

4) Go over Thompkins p. 18 on buoyancy and introduce my [buoyancy notes](https://drive.google.com/file/d/1Id1EECa8qefSD6YhJoJzVp2neEarD9Wa/view?usp=sharing)


#### For Thursday

- Work on {ref}`assign2_hydro_solution`
- Read my buoyancy notes and Thompkins through p. 22

### Lecture 6, week 3, Thursday

1) Continue with buoyancy and Taylor Series using {ref}`worksheet4_taylor_solution`
2) Begin work with cloud models
   - Take a look at 3 cloud model animations
   - Introduce {ref}`xarray_intro`
   - Introduce {ref}`tropical_profiles_solution`

#### For Tuesday

- Read Thompkins through p. 31 on water vapor
- Work on {ref}`assign2_hydro_solution` and {ref}`tropical_profiles_solution`

## Week 4

### Lecture 7, week 4, Tuesday

- From last week:

  - {ref}`worksheet4_taylor_solution`
  - {ref}`assign2_hydro_solution`

- Thompkins treatment of the Clausius-Clapeyron equation in Section 2.1 is pretty terse.  This week
we'll go through it in much more detail, starting with the following material:

  - [Shroeder free energy Chapter 5](https://drive.google.com/file/d/1GrzhzpvWbN0H4dP3Z77M4O88Cyey5k2K/view?usp=sharing)
  - [Maxwell relations notes](https://drive.google.com/file/d/1rJO32NGhTW3TA6noQmSdNl2GNXyrZUgA/view?usp=sharing)
  - [Lohmann Chapter 2 thermodynamics](https://drive.google.com/file/d/1YrHBkCiFduYqjd_Dn8_3fZJocmqva1PG/view?usp=sharing)
  - {ref}`worksheet5_latent_solution`

#### For Thursday

- Read the [Schroeder excerpt](https://drive.google.com/file/d/1GrzhzpvWbN0H4dP3Z77M4O88Cyey5k2K/view?usp=sharing) on
  free energies.
- Read [Lohmann section 2.4](https://drive.google.com/file/d/1YrHBkCiFduYqjd_Dn8_3fZJocmqva1PG/view?usp=sharing) up to p. 48 on the Clausius-Clapeyron equation.

- Work on {ref}`tropical_profiles_solution`.

### Lecture 8, week 4, Thursday

- From last time: {ref}`worksheet5_latent_solution`
- This week's worksheet:  {ref}`worksheet6_gibbs_solution`
- The first part of Assigment 4, due Friday midnight: {ref}`assign4_rootfind_solution`

New material covering Thompkins derivation of the Clausius-Clapeyron equation on page 27

- [Clausius-Clapeyron notes](https://drive.google.com/file/d/144fXb0-of3lpwL0xGZYH16yxGJOmvDEh/view?usp=sharing)

A preview of {ref}`assign4_paper_solution`, due the weekend of Feb. 7.  This requires lecture material from next week.

#### For next Tuesday

- Read the [Clausius-Clapeyron notes](https://drive.google.com/file/d/144fXb0-of3lpwL0xGZYH16yxGJOmvDEh/view?usp=sharing)
- Read Thompkins through the end of Chapter 2 (p. 38)
- Read Thompkins Sec 3.2 p. 48-55 on convective stability

Work on {ref}`tropical_profiles_solution` and {ref}`assign4_rootfind_solution`.

## Week 5

### Lecture 9, week 5, Tuesday

- Reading summary:  The most important concepts in the Thompkins reading are:

1) virtual temperature (problem coming on this later): moist air is
   lighter than dry air because $H_2O$ is lighter than $N_2$
   
2) Different measures of liquid water: $\rho_l$, $r_l$, $q_l$

3) isobaric cooling vs. adiabatic expansion on a tephigram (p. 34)

4) The saturated moist adiabat (eqn 2.44. 2.45, p. 35)

5) Equivalent potential temperature, eq.. 2.49 p 35

6) wet bulb temperature eq. 2.53 on p. 37

7) Normand's construction Fig 2.19 p. 38

8) Moist and dry stability, p. 48 Figure 3.17

9) CAPE and $w_{max}$ eq. 3.3 and 3.4 

- For thermodynamic diagram review -- see [Stull Chapter 5](https://www.eoas.ubc.ca/books/Practical_Meteorology/)

- This week:
  - go over: 
    - {ref}`worksheet6_gibbs_solution`
    - {ref}`tropical_profiles_solution`

- Preview three new sets of notes:
  - [Virtual temperature notes](https://drive.google.com/file/d/1EQI86M4I76ARXDV54VOHKQP3J0-kyyDq/view?usp=sharing)
  - [Static energy and moist adiabats](https://drive.google.com/file/d/1BPTvFijtktG4135T1X8DXwI8l_C8xHEt/view?usp=sharing)
  - [Carnot cycle review](https://drive.google.com/file/d/131QoLZPjajoCtreV6muYOxDlF6jvRH4s/view?usp=sharing)
  
- Worksheet demo (no handin): {ref}`mixing_line_example`
- A [blank tephigram](https://docs.google.com/presentation/d/1tgA2Ia1sWPzZVlgspuW7MOUfokDDZmfvS2eDTtFr-Mw/edit?usp=sharing) to use in {ref}`assign4_paper_solution`

#### For Tuesday

- Read the three new sets of notes and work on {ref}`assign4_rootfind_solution` and {ref}`assign4_paper_solution`

### Lecture 10, week 5, Thursday

- Final assignment before break (note new due dates): {ref}`thompkins_solutions`
- New worksheets:
  - {ref}`convec_lifting`
  - {ref}`moist_adiabats`
  - {ref}`tropical_engine_solution`

#### For next Tuesday 
   - Finish Thompkins reading
   - Work on Assigments 4b and 5


## Week 6

### Lecture 11, week 6, Tuesday

- Go over the {ref}`assign4_rootfind_solution` -- demonstrating list comprehensions in Python
- Here is an initial version of the  [mid-term equation sheet](https://drive.google.com/file/d/1L8am0jnhBdlTNzQNR8lE-Iyw0V6lwqg2/view?usp=sharing)
- Go over Thompkins p. 37 on the wet bulb temperature using my [wet bulb notes](https://drive.google.com/file/d/1wsfzaorhVKiMwDlrolf2gtTrSLT4p4iM/view?usp=sharing)
- The {ref}`cape_part1` notebook shows how to calculate CAPE from a sounding as described in Thompkins p. 49
- Today's worksheet on Taylor series {ref}`worksheet7_taylor_solution`

#### For next Thursday

- read Thompkins pp. 55-60 on mixing and downdrafts
- Work on Assigments 4b and 5

### Lecture 12, week 6, Thursday

- {ref}`worksheet7_taylor_solution`
- [Wet bulb vs. $\theta_e$](https://drive.google.com/file/d/1j8vhmRtGaOlpe0JAAS-fvZmH9H8U2ay2/view?usp=sharing)
- {ref}`mixing_line_calc`
- [2016 midterm](https://drive.google.com/file/d/1Nj31FAXyWur9qffnSflhx6CcLBta4PCj/view?usp=sharing)
- [2012 midterm](https://drive.google.com/file/d/1jFKalQ9MXAvNu824LKVyHtNEleZhgI0q/view?usp=sharing)

## Week 7

- midterm break

## Week 8

### Lecture 17, week 8, Tuesday

- go over {ref}`assign4_paper_solution`
- go over {ref}`thompkins_solutions`
- go over the {ref}`study_topics` for the midterm
- new: {ref}`gibbs_notes`

### Lecture 18, week 8, Thursday

- midterm exam


#### Do for next Tuesday

- review Thompkins pp. 45-63 on convection and mixing,
  especially pp. 56 to 61 on [entrainment and mixing]
  
- read my [notes on entrainment](https://drive.google.com/file/d/1VgY5LWMSmavpoHFQYFE-wkPxGnD5x44O/view?usp=sharing)

## Week 9

### Lecure 19, week 9 Tuesday

- Modeling an entraining plume

- Go over my [notes on entrainment](https://drive.google.com/file/d/1VgY5LWMSmavpoHFQYFE-wkPxGnD5x44O/view?usp=sharing)
and Thompkins pp. 56 to 61

- Introduce {ref}`entrain_cloud`
- Work on {ref}`worksheet8_interp1d_solution`  which is the first part of Assignment 6

#### Do for next Thursday

- Read Thompkins pp. 65-73 on the Köhler equation

### Lecture 20, week 9, Thursday

- {ref}`2026_midterm_solution`

- {ref}`assign6_koehler`  (due next Friday March 13 midnight)

- [Koehler notes](https://drive.google.com/file/d/1SU64h2sGfMsOygqHL1f6tK_K2HgLv1Nr/view?usp=sharing)

#### Do for next Tuesday

- Read [Lohmann Chapter 6.1-6.5](https://drive.google.com/file/d/1vXbVXKs4E4-a-CiBQp-x6gmA_uxJH8I6/view?usp=sharing)
and my entrainment notes and work on Assignment 6

## Week 10

### Lecture 21, week 10, Tuesday

- Review {ref}`gibbs_notes`

- In-class: derive 

  $$
    dg= -\phi dT + v de
  $$
  
  Starting from the first and second law

- Go over

  - Thompkins pp. 56 to 61
  - [Lohmann Chapter 6.1-6.5](https://drive.google.com/file/d/1vXbVXKs4E4-a-CiBQp-x6gmA_uxJH8I6/view?usp=sharing)
  - [Koehler notes](https://drive.google.com/file/d/1SU64h2sGfMsOygqHL1f6tK_K2HgLv1Nr/view?usp=sharing)

- Work on {ref}`worksheet9_koehler_solution` and {ref}`assign6_koehler_solution`

#### For Thursday

Read Thompkins pages 73-75 on condensation growth


### Lecture 22, week 10, Thursday

- Review Thompkins on the droplet growth equation, and my [droplet growth notes](https://drive.google.com/file/d/1mF4ES2IewpnLB-lQHOBuK5JwSnelqkYs/view?usp=sharing)
- How quickly does the vapor concentration adjust to changes, compared to cloud droplets?
  Go over the calculation in my  [timescale notes](https://drive.google.com/file/d/1i_2_LTQ496FScNK8V-b3pb-rvh6Qzapa/view?usp=sharing)
- Start on [Lohmann Aerosols Chapter 5](https://drive.google.com/file/d/1frWt9QyWx8RhSQ8FO3Q6IM4vH4NgrOsU/view) and my
  notes on [aerosol size distributions](https://drive.google.com/file/d/1YYr2dFO4csFNsdOL7IEZJCPEh6LXMoof/view?usp=sharing)

#### For Tuesday

- Read
  - Lohmann Chapter 5 pp. 115-131 on aerosol size distributions
  - Lohmann Chapter 6 pp. 177-182 on Cloud condensation nuclei counters
  - My notes on droplet growth, timescales, and size distributions
  
- Start on {ref}`assign7_aerosols`


## Week 11

### Lecture 23, week 11, Tuesday

- [Final exam practice](https://drive.google.com/file/d/1o9HOkiaXsestzA7P0Njt9X1r5IlEj34D/view?usp=sharing)

- Go over Assigment 6
  - {ref}`worksheet8_interp1d_solution`
  - {ref}`assign6_cape_solution`
  - {ref}`assign6_koehler_solution`
  
- {ref}`worksheet10_aerosol_dists_solution`

#### For Thursday

- Read Thompkins 76-81 on droplets forming precipitation

### Lecture 24, week 11, Thursday

- Review {ref}`worksheet10_aerosol_dists_solution`
- Flashback: [Accurate enthalpy and $\theta_e$](https://drive.google.com/file/d/1kGjCV3tVA4umfoU1KFlkCAbHtXqdOStt/view?usp=sharing)
- Introduce {ref}`assignment_precip` which is based on 
  [Wallace and Hobbs Section 6.4](https://drive.google.com/file/d/1wcG9cKgb6MT7XklFuEZJtcAIQ1uH9mAP/view?usp=sharing)
- {ref}`worksheet11_dropgrow`


#### For Tuesday

Read [Wallace and Hobbs Section 6.4](https://www.dropbox.com/scl/fi/r7zhctw6usest2rirpe2f/wallace_hobbs_chap6.pdf?rlkey=oxdanppj9wkocmwe4nn4cc4fi&dl=0) and work on Assignment 8


## Week 12

### Lecture 24, week 12, Tuesday



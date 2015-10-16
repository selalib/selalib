#--------------------------------------------------------
# Python script that takes a string, normally the name
# of a SeLaLib module, and gives a new string that
# respects the library conventions, typically in the
# form "sll_m_<module_name>".
# To execute : renaming.py (shows some examples)
# or import it and call "renaming(<string>)"
# Author: Laura S. Mendoza
#--------------------------------------------------------

# Globals variables:
# Types of entities to be renamed. Curently script only working with modules.
# In the future there could be additional variables for other entities
# (ie. simulations, functions, subroutines, ...)
RENAME_MOD = 1


def keywords(type = RENAME_MOD):
    """ Returns an array of keywords, for the type of entity specified.
    By default, returns SeLaLib's modules keywords.

    Args:
       type (int): type of entity to be renamed."""

    if (type == RENAME_MOD):
        list_key = []
        list_key+= ["sll"]
        list_key+= ["m"]
        list_key+= ["module"]
        list_key+= ["mod"]
        list_key+= ["sim"]
        list_key+= ["simulation"]
    else:
        raise SystemExit("Error in keywords(). Undefined parameter type="+type)

    return list_key


def convention(splitted, type = RENAME_MOD):
    """ Returns a string common to the type of entity passed in parameter.
    By default, returns SeLaLib's modules string found at the begining.

    Args:
       type (int): type of entity to be renamed."""

    if (type == RENAME_MOD):
        conv = "sll_m"
        for word in splitted :
            if (word == "simulation") or (word == "sim") :
                conv = "sll_m_sim"
    else:
        raise SystemExit("Error in keywords(). Undefined parameter type="+type)

    return conv


def renaming(original_name):
    """ Changes a string to make it respects SeLaLib naming convention

    Args :
       original_name (str) : original module name. """

    list_key   = keywords()
    splitted   = original_name.split("_")
    modulename = convention(splitted)

    for word in splitted :
        if not word in list_key :
            modulename += "_"
            modulename += word

    return modulename


def main():
    print("sll_module_modname      -->    " + renaming("sll_module_modname"))
    print("modname                 -->    " + renaming("modname"))
    print("modname_module          -->    " + renaming("modname_module"))
    print("sll_m_modname           -->    " + renaming("sll_m_modname"))
    print("specific_modname_1d_mod -->    " + renaming("specific_modname_1d_mod"))
    print("simulation_modname_mod  -->    " + renaming("simulation_modname_mod"))


if __name__ == "__main__":
        main()

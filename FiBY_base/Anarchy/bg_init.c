#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "allvars.h"
#include "bg_vars.h"
#include "bg_proto.h"
#include "proto.h"


#if defined(BG_COOLING) || defined(BG_STELLAR_EVOLUTION)

char ElementNames[BG_NELEMENTS][EL_NAME_LENGTH];

double SolarMetallicity, CalciumOverSilicon, SulphurOverSilicon;

void InitChemistry(void)
{
  int Silicon_SPH_Index, Calcium_SPH_Index, Sulphur_SPH_Index;

  /* copy element names to ElementNames */
  strcpy(ElementNames[0], "Hydrogen");
  strcpy(ElementNames[1], "Helium");
  strcpy(ElementNames[2], "Carbon");
  strcpy(ElementNames[3], "Nitrogen");
  strcpy(ElementNames[4], "Oxygen");
  strcpy(ElementNames[5], "Neon");
  strcpy(ElementNames[6], "Magnesium");
  strcpy(ElementNames[7], "Silicon");
  strcpy(ElementNames[8], "Iron");

  /* check we track silicon if calcium or sulphur are scaled with silicon */
  Silicon_SPH_Index = element_index("Silicon");
  Calcium_SPH_Index = element_index("Calcium");
  Sulphur_SPH_Index = element_index("Sulphur");
}

#endif


#if defined(BG_COOLING) || defined(BG_STELLAR_EVOLUTION)

/*
 * ----------------------------------------------------------------------
 * This routine checks whether a given element is tracked
 * ----------------------------------------------------------------------
 */

int element_present(char *element_name)
{
  int i;

  for(i = 0; i < BG_NELEMENTS; i++)
    if(strcmp(ElementNames[i], element_name) == 0)
      return 0;

  /* element not found */
  return -1;
}


/*
 * ----------------------------------------------------------------------
 * This routine returns the index of a given tracked element
 * ----------------------------------------------------------------------
 */

int element_index(char *element_name)
{
  int i;

  for(i = 0; i < BG_NELEMENTS; i++)
    if(strcmp(ElementNames[i], element_name) == 0)
      return i;

  /* element not found */
  return -1;
}


/*
 * ----------------------------------------------------------------------
 * This routine returns the index of an element (element_name) in a
 * generic table of names (table) of given size (size)
 * DON'T USE IT WITH ElementNames!!!!!!!!!!!!
 * ----------------------------------------------------------------------
 */

int get_element_index(char *table[20], int size, char *element_name)
{
  int i;

  for(i = 0; i < size; i++)
    if(strcmp(table[i], element_name) == 0)
      return i;

  /* element not found */
  return -1;
}

#endif

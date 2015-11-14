#include <stdio.h>
using namespace std;

#include <all_sections.h>

int main() {
  for (int i = 0; i < (sizeof(all_sections) / sizeof(all_sections[0])); i++) {
    printf("%s\n", all_sections[i]->name);
  }
  return 0;
}

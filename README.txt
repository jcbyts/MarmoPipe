## MARMOPIPE
I am trying to make sense of the codebase and get rid of the (many) redundant
copies of code that are already present in the example repository.

We don't need to copy "Import" and "MarmoPipe" everytime we use a function from them!

calling addMarmoPipe from matlab will set the paths appropriately for using import and MarmoPipe

Additionally, I'm going to move the external repositories that are required (e.g., Kilosort) out of
MarmoPipe



# Internal S4 generic exports

This roxygen block exists solely to emit `exportMethods` entries for the
S4 generics created in this file. Without explicit `exportMethods`,
methods registered by the package are visible only within the package
namespace and never reach user-facing calls like `rowSums(adgeMatrix)`
at the R console.

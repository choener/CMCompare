import System.Environment (getArgs)
import Text.Printf
import System.Cmd

main = do
  (runlist:idx':[]) <- getArgs
  let idx = (read idx' :: Int) -1 -- lists start at 0
  file <- readFile runlist
  let job = read (lines file !! idx) :: [(String,String)]
  -- execute system call
  let exec (i,j) = printf "./CMCompare testdata/%s testdata/%s -v" i j
  mapM_ (system . exec) job

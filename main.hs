{-# LANGUAGE TypeApplications #-}
import Numeric.LinearAlgebra hiding (Matrix)  -- Импортируем всё, кроме Matrix
import qualified Numeric.LinearAlgebra as LA  -- Используем квалификатор LA
import System.IO  -- Для работы с файлами
import Text.Printf (printf)  -- Для форматированного вывода чисел
import Control.Monad (when)  -- Для проверки условий

-- Определение типа для матриц
type Matrix = LA.Matrix Double  -- Используем квалификатор LA для Matrix из библиотеки

-- Функция для решения уравнения Риккати методом Настрема
solveRiccati :: Matrix -> Matrix -> Matrix -> Double -> Double -> Int -> Int -> Matrix
solveRiccati a b q r dt order iterations = p
  where
    -- Начальное приближение для P(t)
    p0 = LA.konst 0 (rows a, cols a)  -- Создаём нулевую матрицу
    
    -- Функция для вычисления производной P(t)
    dPdt :: Matrix -> Matrix
    dPdt p = aT LA.<> p + p LA.<> a - LA.scale rInv (p LA.<> b LA.<> bT LA.<> p) + q
      where
        aT = tr a
        bT = tr b
        rInv = 1 / r  -- Обратный скаляр
    
    -- Метод Настрема 2-го порядка
    nystrom2 :: Matrix -> Matrix
    nystrom2 p = p + LA.scale dt k2
      where
        k1 = dPdt p
        k2 = dPdt (p + LA.scale (dt/2) k1)
    
    -- Метод Настрема 3-го порядка
    nystrom3 :: Matrix -> Matrix
    nystrom3 p = p + LA.scale (dt/6) (k1 + 4*k2 + k3)
      where
        k1 = dPdt p
        k2 = dPdt (p + LA.scale (dt/2) k1)
        k3 = dPdt (p + LA.scale dt (-k1 + 2*k2))
    
    -- Метод Настрема 4-го порядка
    nystrom4 :: Matrix -> Matrix
    nystrom4 p = p + LA.scale (dt/6) (k1 + 2*k2 + 2*k3 + k4)
      where
        k1 = dPdt p
        k2 = dPdt (p + LA.scale (dt/2) k1)
        k3 = dPdt (p + LA.scale (dt/2) k2)
        k4 = dPdt (p + LA.scale dt k3)
    
    -- Выбор метода в зависимости от порядка
    nystromMethod = case order of
      2 -> nystrom2
      3 -> nystrom3
      4 -> nystrom4
      _ -> error "Поддерживаются только порядки 2, 3 и 4."
    
    -- Итеративное решение
    p = iterate nystromMethod p0 !! iterations

-- Функция для чтения матрицы из файла
readMatrixFromFile :: FilePath -> IO Matrix
readMatrixFromFile path = do
  contents <- readFile path  -- Чтение содержимого файла
  let rows = map (map read . words) (lines contents)  -- Парсинг строк в числа
  return (LA.fromLists rows)  -- Преобразование в матрицу

-- Функция для записи матрицы в файл
writeMatrixToFile :: FilePath -> Matrix -> IO ()
writeMatrixToFile path matrix = do
  let rows = LA.toLists matrix
  let formattedRows = map (unwords . map (printf "%.6f")) rows
  writeFile path (unlines formattedRows)

-- Основная функция
main :: IO ()
main = do
  -- Ввод параметров
  putStrLn "Введите шаг интегрирования (dt):"
  dt <- readLn  -- Чтение шага интегрирования
  
  putStrLn "Введите порядок метода (2, 3 или 4):"
  order <- readLn  -- Чтение порядка метода
  
  putStrLn "Введите количество итераций:"
  iterations <- readLn  -- Чтение количества итераций
  
  -- Проверка корректности ввода
  when (order < 2 || order > 4) $ do
    putStrLn "Ошибка: порядок метода должен быть 2, 3 или 4."
    return ()
  
  when (dt <= 0) $ do
    putStrLn "Ошибка: шаг интегрирования должен быть положительным."
    return ()
  
  when (iterations <= 0) $ do
    putStrLn "Ошибка: количество итераций должно быть положительным."
    return ()
  
  -- Чтение матриц из файлов
  a <- readMatrixFromFile "a.txt"
  b <- readMatrixFromFile "b.txt"
  q <- readMatrixFromFile "q.txt"
  r <- readFile "r.txt" >>= return . read  -- Чтение скаляра r
  
  -- Решение уравнения Риккати
  let p = solveRiccati a b q r dt order iterations
  
  -- Вывод размерности матрицы
  putStrLn $ "Размерность матрицы P: " ++ show (rows p) ++ "x" ++ show (cols p)
  
  -- Сохранение всей матрицы в файл
  writeMatrixToFile "output_matrix.txt" p
  putStrLn "Матрица P сохранена в файл output_matrix.txt"
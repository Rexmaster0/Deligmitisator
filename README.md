Автоматический решатель РК2 по прикладной механике

Для запуска программы необходимо заполнить 3 массива в начале программы.<br/>
В args заполняются коэффициенты перед EI и l для первого и второго конечного элемента.<br/>
В limits заполняется строка с названиями типов ограничений на узлах.<br/>
В forces заполняются силы дейтсвующие на узлы. Не забудьте про знаки у сил.<br/>
![ca0198d3-5451-4999-a6d4-88acbd5fff0d](https://github.com/user-attachments/assets/e522cc68-e7e5-44b3-9fbf-0dd985183edd)
Для этого примера массивы заполняются так:<br/>
args = [4, a, 4, a]<br/>
limits = "поршень, шарнир односвязный, шарнир двусвязный"<br/>
forces = [2*F, 0, 0, 0, 0, 3*F*l]<br/>

После запуска программы файл rk2.tex переформатируется в rk2.pdf, в нем находится все решение кроме рисунков самих стержней.<br/>
Но если у вас не Windows система, то для переноса rk2.tex в pdf вам надо закинуть его на этот сайт:<br/>
https://products.groupdocs.app/conversion/tex-to-pdf"
